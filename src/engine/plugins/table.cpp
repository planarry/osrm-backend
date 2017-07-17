#include "engine/plugins/table.hpp"

#include "engine/api/table_api.hpp"
#include "engine/api/table_parameters.hpp"
#include "engine/routing_algorithms/many_to_many.hpp"
#include "engine/search_engine_data.hpp"
#include "util/json_container.hpp"
#include "util/string_util.hpp"
#include "engine/guidance/leg_geometry.hpp"
#include "engine/guidance/assemble_geometry.hpp"

#include <cstdlib>

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include <boost/assert.hpp>

namespace osrm
{
namespace engine
{
namespace plugins
{

TablePlugin::TablePlugin(datafacade::BaseDataFacade &facade, const int max_locations_distance_table)
    : BasePlugin{facade}, distance_table(&facade, heaps), shortest_path(&facade, heaps),
      max_locations_distance_table(max_locations_distance_table)
{
}

Status TablePlugin::HandleRequest(const api::TableParameters &params, util::json::Object &result)
{
    BOOST_ASSERT(params.IsValid());

    if (!CheckAllCoordinates(params.coordinates)) {
        return Error("InvalidOptions", "Coordinates are invalid", result);
    }

    if (params.bearings.size() > 0 && params.coordinates.size() != params.bearings.size()) {
        return Error(
                "InvalidOptions", "Number of bearings does not match number of coordinates", result);
    }

    // Empty sources or destinations means the user wants all of them included, respectively
    // The ManyToMany routing algorithm we dispatch to below already handles this perfectly.
    const auto num_sources =
            params.sources.empty() ? params.coordinates.size() : params.sources.size();
    const auto num_destinations =
            params.destinations.empty() ? params.coordinates.size() : params.destinations.size();

    if (max_locations_distance_table > 0 &&
        ((num_sources * num_destinations) >
         static_cast<std::size_t>(max_locations_distance_table * max_locations_distance_table))) {
        return Error("TooBig", "Too many table coordinates", result);
    }

    auto snapped_phantoms = SnapPhantomNodes(GetPhantomNodes(params));
    auto result_table = distance_table(snapped_phantoms, params.sources, params.destinations, params.graph_flag);

    if (result_table.first.first.empty()) {
        return Error("NoTable", "No table found", result);
    }

    api::TableAPI table_api{facade, params};
    table_api.MakeResponse(result_table.first, snapped_phantoms, result);

    std::pair<unsigned, unsigned> time_limits = facade.GetLimitsOfTime();

    if ((params.time_period_from) && (params.time_period_to > params.time_period_from)
        && (time_limits.first < params.time_period_to)
        && (time_limits.second > params.time_period_from)) {
        std::vector<PhantomNodes> nodes;

        std::vector<unsigned long> our_sources = params.sources, our_destinations = params.destinations;
        if (our_sources.empty()) {
            our_sources.reserve(snapped_phantoms.size());
            for (unsigned long i = 0; i < snapped_phantoms.size(); ++i) {
                our_sources.push_back(i);
            }
        }
        if (our_destinations.empty()) {
            our_destinations.reserve(snapped_phantoms.size());
            for (unsigned long i = 0; i < snapped_phantoms.size(); ++i) {
                our_destinations.push_back(i);
            }
        }

        // Create result table
        std::vector<std::vector<std::vector<api::TableAPI::AdditionalWeights>>> addition_table_result;
        addition_table_result.resize(our_sources.size());
        for (unsigned k = 0; k < our_sources.size(); ++k) {
            addition_table_result[k].resize(our_destinations.size());
        }

        for (size_t index_source = 0; index_source < our_sources.size(); ++index_source) {
            for (size_t index_destination = 0; index_destination < our_destinations.size(); ++index_destination) {
                api::TableAPI::AdditionalWeights current_result{params.time_period_from, params.time_period_to, 0};
                // Check if it is the same node
                if (our_sources[index_source] == our_destinations[index_destination]) {
                    addition_table_result[index_source][index_destination].push_back(current_result);
                    continue;
                }
                // Check if both nodes are on same edge
                unsigned part_of_edge_weight = 0;
                if (snapped_phantoms[our_sources[index_source]].forward_segment_id.id ==
                    snapped_phantoms[our_destinations[index_destination]].forward_segment_id.id
                    && snapped_phantoms[our_sources[index_source]].reverse_segment_id.id ==
                       snapped_phantoms[our_destinations[index_destination]].reverse_segment_id.id) // one way
                {
                    part_of_edge_weight = (unsigned) abs(
                            snapped_phantoms[our_sources[index_source]].GetForwardWeightPlusOffset()
                            - snapped_phantoms[our_destinations[index_destination]].GetForwardWeightPlusOffset());
                } else if (snapped_phantoms[our_sources[index_source]].reverse_segment_id.id ==
                           snapped_phantoms[our_destinations[index_destination]].forward_segment_id.id
                           && snapped_phantoms[our_sources[index_source]].forward_segment_id.id ==
                              snapped_phantoms[our_destinations[index_destination]].reverse_segment_id.id) // opposite way
                {
                    part_of_edge_weight = (unsigned) abs(
                            snapped_phantoms[our_sources[index_source]].GetForwardWeightPlusOffset()
                            - snapped_phantoms[our_destinations[index_destination]].GetReverseWeightPlusOffset());
                }

                nodes.clear();
                nodes.push_back({snapped_phantoms[our_sources[index_source]],
                                 snapped_phantoms[our_destinations[index_destination]]});

                InternalRouteResult raw_route;
                shortest_path(nodes, boost::optional<bool>(true), raw_route);

                if (raw_route.is_valid()) {
                    osrm::engine::guidance::LegGeometry legGeometry = osrm::engine::guidance::assembleGeometry(
                            facade,
                            raw_route.unpacked_path_segments.front(),
                            snapped_phantoms[our_sources[index_source]],
                            snapped_phantoms[our_destinations[index_destination]]);

                    // Time cycle
                    unsigned delta_time = 5 * 60;
                    std::vector<unsigned int> started_time;
                    std::vector<unsigned int> result_time;
                    for (unsigned int time_index = params.time_period_from;
                         time_index < params.time_period_to
                                      - result_table.first.first[index_source * our_destinations.size() + index_destination];
                         time_index += delta_time) {
                        started_time.push_back(time_index);
                        result_time.push_back(time_index);
                    }
                    if (started_time.empty())
                        continue;
                    // Correct source and target OSMNodeID
                    //TODO probably we know what node we should put
                    std::pair<OSMNodeID, OSMNodeID> osm_node_start, osm_node_finish;
                    std::vector<NodeID> source_geometry;
                    facade.GetUncompressedGeometry(
                            snapped_phantoms[our_sources[index_source]].reverse_packed_geometry_id, source_geometry);
                    osm_node_start.first = facade.GetOSMNodeIDOfNode(source_geometry[source_geometry.size() -
                                                                                     snapped_phantoms[our_sources[index_source]].fwd_segment_position -
                                                                                     1]);
                    facade.GetUncompressedGeometry(
                            snapped_phantoms[our_sources[index_source]].forward_packed_geometry_id, source_geometry);
                    osm_node_start.second = facade.GetOSMNodeIDOfNode(
                            source_geometry[snapped_phantoms[our_sources[index_source]].fwd_segment_position]);

                    std::vector<NodeID> target_geometry;
                    facade.GetUncompressedGeometry(
                            snapped_phantoms[our_destinations[index_destination]].reverse_packed_geometry_id,
                            target_geometry);
                    osm_node_finish.first = facade.GetOSMNodeIDOfNode(target_geometry[target_geometry.size() -
                                                                                      snapped_phantoms[our_destinations[index_destination]].fwd_segment_position -
                                                                                      1]);
                    facade.GetUncompressedGeometry(
                            snapped_phantoms[our_destinations[index_destination]].forward_packed_geometry_id,
                            target_geometry);
                    osm_node_finish.second = facade.GetOSMNodeIDOfNode(
                            target_geometry[snapped_phantoms[our_destinations[index_destination]].fwd_segment_position]);

                    if (legGeometry.osm_node_ids.size() > 2) {
                        legGeometry.osm_node_ids[0] = (legGeometry.osm_node_ids[1] == osm_node_start.second) ?
                                                      osm_node_start.first : osm_node_start.second;
                        legGeometry.osm_node_ids[legGeometry.osm_node_ids.size() - 1] =
                                (legGeometry.osm_node_ids[legGeometry.osm_node_ids.size() - 2] ==
                                 osm_node_finish.second) ?
                                osm_node_finish.first : osm_node_finish.second;
                    }

                    // Check every edge (a pair of nodes)
                    if ((legGeometry.annotations.size() ==
                         legGeometry.osm_node_ids.size() - 1)  // each edge(pair of nodes) has a weight
                        && !part_of_edge_weight) {
                        for (unsigned curr_node = 0; curr_node < legGeometry.annotations.size(); ++curr_node) {
                            std::vector<SegmentAddition>::iterator addition_time_from;
                            std::vector<SegmentAddition>::iterator addition_time_to;
                            bool has_addition_time = facade.GetIteratorsOfAdditionWeights(
                                    legGeometry.osm_node_ids[curr_node],
                                    legGeometry.osm_node_ids[curr_node + 1],
                                    params.time_period_from,
                                    params.time_period_to,
                                    addition_time_from,
                                    addition_time_to);

                            for (unsigned i = 0; i < started_time.size(); ++i) {
                                bool is_addition = false;
                                if (has_addition_time) {
                                    while ((addition_time_from != addition_time_to)
                                           && (result_time[i] > addition_time_from->time_from)) {
                                        if ((result_time[i] >= addition_time_from->time_from)
                                            && (result_time[i] <= addition_time_from->time_to)) {

                                            is_addition = true;
                                            if (curr_node && curr_node + 1 != legGeometry.annotations.size()) {
                                                if (addition_time_from->speed < legGeometry.annotations[curr_node].duration * 10)
                                                    result_time[i] += legGeometry.annotations[curr_node].duration * 10; // this shouldn't happen at all
                                                else
                                                    result_time[i] += addition_time_from->speed; // actually it is weight
                                            }
                                            else
                                                result_time[i] += addition_time_from->speed
                                                                  * (legGeometry.annotations[curr_node].duration * 10)
                                                                  / (double) (snapped_phantoms[our_sources[index_source]].forward_weight
                                                                              + snapped_phantoms[our_sources[index_source]].reverse_weight);

                                            break;
                                        }
                                        addition_time_from++;
                                    }
                                }
                                if (!is_addition)
                                    result_time[i] += (unsigned int) (legGeometry.annotations[curr_node].duration * 10);
                                if (i && (result_time[i] < result_time[i - 1])) // element i can't outdrive i-1
                                    result_time[i] = result_time[i - 1];
                            }
                        }
                    }

                    // Convert results to JSON structure
                    if (!part_of_edge_weight) {
                        // Take first element
                        result_time[0] -= started_time[0];
                        current_result = (api::TableAPI::AdditionalWeights) {started_time[0],
                                                                             started_time[0] + (delta_time - 1),
                                                                             result_time[0]};
                        for (unsigned i = 1; i < started_time.size(); ++i) {
                            result_time[i] -= started_time[i];
                            if (current_result.additional_weight ==
                                result_time[i]) { // TODO better to check using correction like r>=w && r<w+correction
                                current_result.start_time_to = started_time[i] + (delta_time - 1);
                            } else {
                                addition_table_result[index_source][index_destination].push_back(current_result);
                                current_result.start_time_from = started_time[i];
                                current_result.start_time_to = started_time[i] + (delta_time - 1);
                                current_result.additional_weight = result_time[i];
                            }
                        }
                        addition_table_result[index_source][index_destination].push_back(current_result);
                    } else {
                        current_result.additional_weight = (unsigned) (part_of_edge_weight);
                        addition_table_result[index_source][index_destination].push_back(current_result);
                    }
                }
            }
        }
        table_api.AppendResponse(addition_table_result, result);
    } else { // empty data for math
        std::vector<std::vector<std::vector<api::TableAPI::AdditionalWeights>>> addition_table_result;
        table_api.AppendResponse(addition_table_result, result);
    }
    if (params.graph_flag == 1)
        table_api.AppendGraphResponse(result_table.second, snapped_phantoms, result);

    return Status::Ok;
}
}
}
}
