//
// Created by dan1992 on 31.10.16.
//

#ifndef OSRM_ADDITIONAL_WEIGHT_HPP
#define OSRM_ADDITIONAL_WEIGHT_HPP

#include "typedefs.hpp"

#include <vector>

#include <boost/spirit/include/qi.hpp>
#include <boost/filesystem/fstream.hpp>

#include "util/exception.hpp"
#include "util/simple_logger.hpp"

#include <tbb/concurrent_unordered_map.h>
#include <tbb/spin_mutex.h>
#include <tbb/parallel_for.h>

using namespace osrm;

namespace std
{

    template <> struct hash<std::pair<OSMNodeID, OSMNodeID>>
    {
        std::size_t operator()(const std::pair<OSMNodeID, OSMNodeID> &k) const noexcept
        {
            return static_cast<uint64_t>(k.first) ^ (static_cast<uint64_t>(k.second) << 12);
        }
    };

    template <> struct hash<std::tuple<OSMNodeID, OSMNodeID, OSMNodeID>>
    {
        std::size_t operator()(const std::tuple<OSMNodeID, OSMNodeID, OSMNodeID> &k) const noexcept
        {
            std::size_t seed = 0;
            boost::hash_combine(seed, static_cast<uint64_t>(std::get<0>(k)));
            boost::hash_combine(seed, static_cast<uint64_t>(std::get<1>(k)));
            boost::hash_combine(seed, static_cast<uint64_t>(std::get<2>(k)));
            return seed;
        }
    };
}

// Utilities for LoadEdgeExpandedGraph to restore my sanity
namespace
{
    struct SegmentAddition final
    {
        unsigned time_from = 0, time_to = 0, speed = 0;

        SegmentAddition() {}
        SegmentAddition(unsigned f, unsigned t, unsigned s) : time_from(f), time_to(t), speed(s){}
    };

    struct Segment final
    {
        OSMNodeID from, to;
    };

    struct SpeedSource final
    {
        unsigned speed;
        std::uint8_t source;
        std::vector<SegmentAddition> segments;
    };

    struct SegmentSpeedSource final
    {
        Segment segment;
        SpeedSource speed_source;
    };

    using SegmentSpeedSourceFlatMap = std::vector<SegmentSpeedSource>;

// Binary Search over a flattened key,val Segment storage
    template <typename T = SegmentSpeedSource, template <typename, typename> class Container = std::vector>
    typename Container<T, std::allocator<T>>::iterator find(Container<T, std::allocator<T>> &map, const Segment &key)
    {
        const auto last = map.end();
//        const auto last = end(map);

        const auto by_segment = [](const T &lhs, const T &rhs) {
            return std::tie(lhs.segment.from, lhs.segment.to) >
                   std::tie(rhs.segment.from, rhs.segment.to);
        };

        auto it = std::lower_bound(map.begin(), last, T{key, {0, 0}}, by_segment);
//        auto it = std::lower_bound(begin(map), last, SegmentSpeedSource{key, {0, 0}}, by_segment);

        if (it != last && (std::tie(it->segment.from, it->segment.to) == std::tie(key.from, key.to)))
            return it;

        return last;
    }

// Convenience aliases. TODO: make actual types at some point in time.
// TODO: turn penalties need flat map + binary search optimization, take a look at segment speeds

    using Turn = std::tuple<OSMNodeID, OSMNodeID, OSMNodeID>;
    using TurnHasher = std::hash<Turn>;
    using PenaltySource = std::pair<double, std::uint8_t>;
    using TurnPenaltySourceMap = tbb::concurrent_unordered_map<Turn, PenaltySource, TurnHasher>;

// Functions for parsing files and creating lookup tables

    SegmentSpeedSourceFlatMap
    parse_segment_lookup_from_csv_files(const std::vector<std::string> &segment_speed_filenames)
    {
        // TODO: shares code with turn penalty lookup parse function

        using Mutex = tbb::spin_mutex;

        // Loaded and parsed in parallel, at the end we combine results in a flattened map-ish view
        SegmentSpeedSourceFlatMap flatten;
        Mutex flatten_mutex;

        const auto parse_segment_speed_file = [&](const std::size_t idx) {
            const auto file_id = idx + 1; // starts at one, zero means we assigned the weight
            const auto filename = segment_speed_filenames[idx];

            std::ifstream segment_speed_file{filename, std::ios::binary};
            if (!segment_speed_file)
                throw util::exception{"Unable to open segment speed file " + filename};

            SegmentSpeedSourceFlatMap local;

            std::uint64_t from_node_id{};
            std::uint64_t to_node_id{};
            unsigned speed{};

            for (std::string line; std::getline(segment_speed_file, line);)
            {
                using namespace boost::spirit::qi;

                auto it = begin(line);
                const auto last = end(line);

//            std::string add_data;

                // The ulong_long -> uint64_t will likely break on 32bit platforms
                const auto ok =
                        parse(it,
                              last,                                                                  //
                              (ulong_long >> ',' >> ulong_long >> ',' >> uint_ /*>> *(',' >> *char_)*/), //
                              from_node_id,
                              to_node_id,
                              speed/*,
                      add_data*/); //

                if (!ok /*|| it != last*/)
                    throw util::exception{"Segment speed file " + filename + " malformed"};

                SegmentSpeedSource val{{OSMNodeID{from_node_id}, OSMNodeID{to_node_id}},
                                       {speed, static_cast<std::uint8_t>(file_id)}};

                std::uint64_t from_time{};
                std::uint64_t to_time{};
                unsigned speed_add_time{};

                while(parse(it,
                            last,                                                                  //
                            (',' >> ulong_long >> ',' >> ulong_long >> ',' >> uint_ ), //
                            from_time,
                            to_time,
                            speed_add_time)){
                    if(speed_add_time > 0 && speed_add_time < speed /* * 0.9 */) // probably should use time correction to avoid small deviation
                        val.speed_source.segments.emplace_back(
                                from_time,
                                to_time,
                                speed_add_time);
                }
                local.push_back(std::move(val));
            }

            util::SimpleLogger().Write() << "Loaded speed file " << filename << " with " << local.size()
                                         << " speeds";

            {
                Mutex::scoped_lock _{flatten_mutex};

                flatten.insert(end(flatten),
                               std::make_move_iterator(begin(local)),
                               std::make_move_iterator(end(local)));
            }
        };

        tbb::parallel_for(std::size_t{0}, segment_speed_filenames.size(), parse_segment_speed_file);

        // With flattened map-ish view of all the files, sort and unique them on from,to,source
        // The greater '>' is used here since we want to give files later on higher precedence
        const auto sort_by = [](const SegmentSpeedSource &lhs, const SegmentSpeedSource &rhs) {
            return std::tie(lhs.segment.from, lhs.segment.to, lhs.speed_source.source) >
                   std::tie(rhs.segment.from, rhs.segment.to, rhs.speed_source.source);
        };

        std::stable_sort(begin(flatten), end(flatten), sort_by);

        // Unique only on from,to to take the source precedence into account and remove duplicates
        const auto unique_by = [](const SegmentSpeedSource &lhs, const SegmentSpeedSource &rhs) {
            return std::tie(lhs.segment.from, lhs.segment.to) ==
                   std::tie(rhs.segment.from, rhs.segment.to);
        };

        const auto it = std::unique(begin(flatten), end(flatten), unique_by);

        flatten.erase(it, end(flatten));

        util::SimpleLogger().Write() << "In total loaded " << segment_speed_filenames.size()
                                     << " speed file(s) with a total of " << flatten.size()
                                     << " unique values";

        return flatten;
    }

    void saveAddWeightsToFile(SegmentSpeedSourceFlatMap &flatMap, std::string fileName){
        boost::filesystem::ofstream output_stream(fileName, std::ios::binary);
        size_t countEdges = 0; //+ flatMap.size();
        output_stream.write((char *)&countEdges, sizeof(size_t));
        for (const auto &item : flatMap) {
            size_t countSegments = item.speed_source.segments.size();
            if (countSegments) {
                output_stream.write((char *)&(item.segment), sizeof(Segment));
                output_stream.write((char *)&countSegments, sizeof(size_t));
                output_stream.write((char *)&(item.speed_source.segments[0]), sizeof(SegmentAddition) * countSegments);
                countEdges++;
            }
        }
        output_stream.seekp(std::ios_base::beg);
        output_stream.write((char *)&countEdges, sizeof(size_t));
        output_stream.seekp(std::ios_base::end);
        output_stream.close();
    }

    TurnPenaltySourceMap
    parse_turn_penalty_lookup_from_csv_files(const std::vector<std::string> &turn_penalty_filenames)
    {
        // TODO: shares code with turn penalty lookup parse function
        TurnPenaltySourceMap map;

        const auto parse_turn_penalty_file = [&](const std::size_t idx) {
            const auto file_id = idx + 1; // starts at one, zero means we assigned the weight
            const auto filename = turn_penalty_filenames[idx];

            std::ifstream turn_penalty_file{filename, std::ios::binary};
            if (!turn_penalty_file)
                throw util::exception{"Unable to open turn penalty file " + filename};

            std::uint64_t from_node_id{};
            std::uint64_t via_node_id{};
            std::uint64_t to_node_id{};
            double penalty{};

            for (std::string line; std::getline(turn_penalty_file, line);)
            {
                using namespace boost::spirit::qi;

                auto it = begin(line);
                const auto last = end(line);

                // The ulong_long -> uint64_t will likely break on 32bit platforms
                const auto ok = parse(it,
                                      last, //
                                      (ulong_long >> ',' >> ulong_long >> ',' >> ulong_long >> ',' >>
                                                  double_ >> *(',' >> *char_)), //
                                      from_node_id,
                                      via_node_id,
                                      to_node_id,
                                      penalty); //

                if (!ok || it != last)
                    throw util::exception{"Turn penalty file " + filename + " malformed"};

                map[std::make_tuple(
                        OSMNodeID{from_node_id}, OSMNodeID{via_node_id}, OSMNodeID{to_node_id})] =
                        std::make_pair(penalty, file_id);
            }
        };

        tbb::parallel_for(std::size_t{0}, turn_penalty_filenames.size(), parse_turn_penalty_file);

        return map;
    }
} // anon ns

#endif //OSRM_ADDITIONAL_WEIGHT_HPP
