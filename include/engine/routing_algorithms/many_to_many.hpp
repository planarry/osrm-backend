#ifndef MANY_TO_MANY_ROUTING_HPP
#define MANY_TO_MANY_ROUTING_HPP

#include "engine/routing_algorithms/routing_base.hpp"
#include "engine/search_engine_data.hpp"
#include "util/typedefs.hpp"

#include <boost/assert.hpp>

#include <limits>
#include <memory>
#include <unordered_map>
#include <vector>

namespace osrm
{
namespace engine
{
namespace routing_algorithms
{

template <class DataFacadeT>
class ManyToManyRouting final
    : public BasicRoutingInterface<DataFacadeT, ManyToManyRouting<DataFacadeT>>
{
    using HierarchyTree = std::vector<std::unordered_map<NodeID, NodeID>>;
    using EdgeData = typename DataFacadeT::EdgeData;
    using super = BasicRoutingInterface<DataFacadeT, ManyToManyRouting<DataFacadeT>>;
    using QueryHeapL = SearchEngineData::QueryHeapL;
    SearchEngineData &engine_working_data;

    struct NodeBucket
    {
        unsigned target_id; // essentially a row in the distance matrix
        EdgeWeight distance;
        EdgeLength length;
        NodeBucket(const unsigned target_id, const EdgeWeight distance, const EdgeLength length)
            : target_id(target_id), distance(distance), length(length)
        {
        }
    };

    // FIXME This should be replaced by an std::unordered_multimap, though this needs benchmarking
    using SearchSpaceWithBuckets = std::unordered_map<NodeID, std::vector<NodeBucket>>;

  public:
    ManyToManyRouting(DataFacadeT *facade, SearchEngineData &engine_working_data)
        : super(facade), engine_working_data(engine_working_data)
    {
    }

    std::pair<std::pair<std::vector<EdgeWeight>, std::vector<EdgeLength>>, std::vector<std::list<PointID>>> operator()(
            const std::vector<PhantomNode> &phantom_nodes,
            const std::vector<std::size_t> &source_indices,
            const std::vector<std::size_t> &target_indices,
            const unsigned int needGraph) const
    {
        const unsigned number_of_sources =
            source_indices.empty() ? phantom_nodes.size() : source_indices.size();
        const unsigned number_of_targets =
            target_indices.empty() ? phantom_nodes.size() : target_indices.size();
        const unsigned number_of_entries = number_of_sources * number_of_targets;
        std::pair< std::pair<std::vector<EdgeWeight>, std::vector<EdgeLength>>,
                   std::vector<std::list<PointID>> > result_table(
                std::make_pair(
                std::make_pair(std::vector<EdgeWeight>(number_of_entries, std::numeric_limits<EdgeWeight>::max()),
                               std::vector<EdgeLength>(number_of_entries, std::numeric_limits<EdgeLength>::max())),
                std::vector<std::list<PointID>>(number_of_sources, std::numeric_limits<std::list<PointID>>::max())));
        std::vector<NodeID> cross_nodes_table(number_of_entries, std::numeric_limits<NodeID>::max());

        engine_working_data.InitializeOrClearLengthThreadLocalStorage(
            super::facade->GetNumberOfNodes());

        QueryHeapL &query_heap = *(engine_working_data.forward_heap_l_1);

        SearchSpaceWithBuckets search_space_with_buckets;

        HierarchyTree backward_hierarchy_tree(number_of_targets), forward_hierarchy_tree(number_of_sources);

        std::vector<std::set<NodeID>> phantomes_for_point(number_of_targets);
        std::set<NodeID> all_phantomes;

        unsigned column_idx = 0;
        const auto search_target_phantom = [&](const PhantomNode &phantom) {
            query_heap.Clear();
            // insert target(s) at distance 0

            if (phantom.forward_segment_id.enabled)
            {
                if (needGraph == 1) {
                    phantomes_for_point[column_idx].insert(phantom.forward_segment_id.id);
                    all_phantomes.insert(phantom.forward_segment_id.id);
                }
                query_heap.Insert(phantom.forward_segment_id.id,
                                  std::pair<int, int>(phantom.GetForwardWeightPlusOffset(),
                                                      phantom.GetForwardLengthPlusOffset()),
                                  phantom.forward_segment_id.id);
            }
            if (phantom.reverse_segment_id.enabled)
            {
                if (needGraph == 1) {
                    phantomes_for_point[column_idx].insert(phantom.reverse_segment_id.id);
                    all_phantomes.insert(phantom.reverse_segment_id.id);
                }
                query_heap.Insert(phantom.reverse_segment_id.id,
                                  std::pair<int, int>(phantom.GetReverseWeightPlusOffset(),
                                                      phantom.GetReverseLengthPlusOffset()),
                                  phantom.reverse_segment_id.id);
            }

            // explore search space
            while (!query_heap.Empty())
            {
                BackwardRoutingStep(column_idx, query_heap, search_space_with_buckets, backward_hierarchy_tree, needGraph);
            }
            ++column_idx;
        };

        // for each source do forward search
        unsigned row_idx = 0;
        const auto search_source_phantom = [&](const PhantomNode &phantom) {
            query_heap.Clear();
            // insert target(s) at distance 0

            if (phantom.forward_segment_id.enabled)
            {
                query_heap.Insert(phantom.forward_segment_id.id,
                                  std::pair<int, int>(-phantom.GetForwardWeightPlusOffset(),
                                                      -phantom.GetForwardLengthPlusOffset()),
                                  phantom.forward_segment_id.id);
            }
            if (phantom.reverse_segment_id.enabled)
            {
                query_heap.Insert(phantom.reverse_segment_id.id,
                                  std::pair<int, int>(-phantom.GetReverseWeightPlusOffset(),
                                                      -phantom.GetReverseLengthPlusOffset()),
                                  phantom.reverse_segment_id.id);
            }

            // explore search space
            while (!query_heap.Empty())
            {
                ForwardRoutingStep(row_idx,
                                   number_of_targets,
                                   query_heap,
                                   search_space_with_buckets,
                                   result_table,
                                   forward_hierarchy_tree,
                                   cross_nodes_table,
                                   needGraph);
            }
            ++row_idx;
        };

        if (target_indices.empty())
        {
            for (const auto &phantom : phantom_nodes)
            {
                search_target_phantom(phantom);
            }
        }
        else
        {
            for (const auto index : target_indices)
            {
                const auto &phantom = phantom_nodes[index];
                search_target_phantom(phantom);
            }
        }

        if (source_indices.empty())
        {
            for (const auto &phantom : phantom_nodes)
            {
                search_source_phantom(phantom);
            }
        }
        else
        {
            for (const auto index : source_indices)
            {
                const auto &phantom = phantom_nodes[index];
                search_source_phantom(phantom);
            }
        }

        if (needGraph == 1)
            BuildFullGraph(result_table,
                           cross_nodes_table,
                           forward_hierarchy_tree,
                           backward_hierarchy_tree,
                           phantomes_for_point,
                           all_phantomes,
                           number_of_sources,
                           number_of_targets);

        return result_table;
    }

    void ForwardRoutingStep(const unsigned row_idx,
                            const unsigned number_of_targets,
                            QueryHeapL &query_heap,
                            const SearchSpaceWithBuckets &search_space_with_buckets,
                            std::pair<std::pair<std::vector<EdgeWeight>, std::vector<EdgeLength>>,
                                    std::vector<std::list<PointID>>> &result_table,
                            HierarchyTree &hierarchy_tree, std::vector<NodeID> &cross_nodes_table,
                            const unsigned int needGraph) const
    {
        const NodeID node = query_heap.DeleteMin();
        const int source_distance = query_heap.GetKey(node).first;
        const int source_length = query_heap.GetKey(node).second;

        if (needGraph == 1)
            hierarchy_tree[row_idx][node] = query_heap.GetData(node).parent;
        // check if each encountered node has an entry
        const auto bucket_iterator = search_space_with_buckets.find(node);
        // iterate bucket if there exists one
        if (bucket_iterator != search_space_with_buckets.end())
        {
            const std::vector<NodeBucket> &bucket_list = bucket_iterator->second;
            for (const NodeBucket &current_bucket : bucket_list)
            {
                // get target id from bucket entry
                const unsigned column_idx = current_bucket.target_id;
                const int target_distance = current_bucket.distance;
                const int target_length = current_bucket.length;
                auto &current_distance = result_table.first.first[row_idx * number_of_targets + column_idx];
                auto &current_length = result_table.first.second[row_idx * number_of_targets + column_idx];
                // check if new distance is better
                const EdgeWeight new_distance = source_distance + target_distance;
                const int new_length = source_length + target_length;
                if (new_distance < 0)
                {
                    const EdgeWeight loop_weight = super::GetLoopWeight(node);
                    const int new_distance_with_loop = new_distance + loop_weight;
                    const int loop_length = super::GetLoopLength(node);
                    const int new_length_with_loop = new_length + loop_length;
                    if (loop_weight != INVALID_EDGE_WEIGHT && new_distance_with_loop >= 0)
                    {
                        current_distance = std::min(current_distance, new_distance_with_loop);
                        current_length = std::min(current_length, new_length_with_loop);
                    }
                }
                else if (new_distance < current_distance)
                {
                    if (needGraph == 1)
                        cross_nodes_table[row_idx * number_of_targets + column_idx] = node;
                    result_table.first.first[row_idx * number_of_targets + column_idx] = new_distance;
                    result_table.first.second[row_idx * number_of_targets + column_idx] = new_length;
                }
            }
        }
        if (StallAtNode<true>(node, source_distance, query_heap))
        {
            return;
        }
        RelaxOutgoingEdges<true>(node, source_distance, source_length, query_heap);
    }

    void BackwardRoutingStep(const unsigned column_idx,
                             QueryHeapL &query_heap,
                             SearchSpaceWithBuckets &search_space_with_buckets,
                             HierarchyTree &hierarchy_tree,
                             const unsigned int needGraph) const
    {
        const NodeID node = query_heap.DeleteMin();
        const int target_distance = query_heap.GetKey(node).first;
        const int target_length = query_heap.GetKey(node).second;

        if (needGraph == 1)
            hierarchy_tree[column_idx][node] = query_heap.GetData(node).parent;
        // store settled nodes in search space bucket
        search_space_with_buckets[node].emplace_back(column_idx, target_distance, target_length);

        if (StallAtNode<false>(node, target_distance, query_heap))
        {
            return;
        }

        RelaxOutgoingEdges<false>(node, target_distance, target_length, query_heap);
    }

    template <bool forward_direction>
    inline void
    RelaxOutgoingEdges(const NodeID node, const EdgeWeight distance, const int length,
                       QueryHeapL &query_heap) const
    {
        for (auto edge : super::facade->GetAdjacentEdgeRange(node))
        {
            const auto &data = super::facade->GetEdgeData(edge);
            const bool direction_flag = (forward_direction ? data.forward : data.backward);
            if (direction_flag)
            {
                const NodeID to = super::facade->GetTarget(edge);
                const int edge_weight = data.distance;
                const int edge_length = (int) data.length;

                BOOST_ASSERT_MSG(edge_weight > 0, "edge_weight invalid");
                const int to_distance = distance + edge_weight;
                const int to_length = length + edge_length;

                // New Node discovered -> Add to Heap + Node Info Storage
                if (!query_heap.WasInserted(to))
                {
                    query_heap.Insert(to, std::pair<int, int>(to_distance, to_length), node);
                }
                // Found a shorter Path -> Update distance
                else if (to_distance < query_heap.GetKey(to).first)
                {
                    // new parent
                    query_heap.GetData(to).parent = node;
                    query_heap.DecreaseKey(to, std::pair<int, int>(to_distance, to_length));
                }
            }
        }
    }

    // Stalling
    template <bool forward_direction>
    inline bool
    StallAtNode(const NodeID node, const EdgeWeight distance, QueryHeapL &query_heap) const
    {
        for (auto edge : super::facade->GetAdjacentEdgeRange(node))
        {
            const auto &data = super::facade->GetEdgeData(edge);
            const bool reverse_flag = ((!forward_direction) ? data.forward : data.backward);
            if (reverse_flag)
            {
                const NodeID to = super::facade->GetTarget(edge);
                const int edge_weight = data.distance;
                BOOST_ASSERT_MSG(edge_weight > 0, "edge_weight invalid");
                if (query_heap.WasInserted(to))
                {
                    if (query_heap.GetKey(to).first + edge_weight < distance)
                    {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    void BuildFullGraph(std::pair<std::pair<std::vector<EdgeWeight>, std::vector<EdgeLength>>,
                                std::vector<std::list<PointID>>> &result_table,
                        const std::vector<NodeID> &cross_nodes_table,
                        const HierarchyTree &forward_hierarchy_tree,
                        const HierarchyTree &backward_hierarchy_tree,
                        const std::vector<std::set<NodeID>> &phantomes_for_point,
                        const std::set<NodeID> &all_phantomes,
                        const unsigned number_of_sources,
                        const unsigned number_of_targets) const {
        std::map<std::pair<NodeID, NodeID>, std::set<NodeID>> unpack_cache;
        for (PointID i = 0; i < number_of_sources; ++i) {
            std::vector<PointID> idxs;
            for (PointID j = 0; j < number_of_targets; ++j)
                if (i != j) idxs.push_back(j);
            std::sort(idxs.begin(), idxs.end(), [i, number_of_targets, &result_table](PointID l, PointID r){
                return result_table.first.first[i * number_of_targets + l] < result_table.first.first[i * number_of_targets + r];
            });

            std::set<NodeID> blocked_cross_nodes, all_blocked, reached_target_phantomes;
            for (PointID j = 0; j < std::min(number_of_targets - 1, 100u); ++j) {
                NodeID cur_node = cross_nodes_table[i * number_of_targets + idxs[j]];
                if (blocked_cross_nodes.find(cur_node) != blocked_cross_nodes.end()) {
//                    SimpleLogger().Write()<<"Bloking "<<i<<" "<<idxs[j]<<". blocked_cross_nodes.";
                    continue;
                }

                all_blocked.clear();
                std::set<NodeID> temp_set;
                std::set_union(phantomes_for_point[i].begin(), phantomes_for_point[i].end(),
                               phantomes_for_point[idxs[j]].begin(), phantomes_for_point[idxs[j]].end(),
                               std::inserter(temp_set, temp_set.begin()));
                std::set_difference(all_phantomes.begin(), all_phantomes.end(),
                                    temp_set.begin(), temp_set.end(),
                                    std::inserter(all_blocked, all_blocked.begin()));
                all_blocked.insert(blocked_cross_nodes.begin(), blocked_cross_nodes.end());

                bool is_first_node = true;
                while(forward_hierarchy_tree[i].count(cur_node) && cur_node != forward_hierarchy_tree[i].at(cur_node))
                {
                    if(SetHasIntersection(all_blocked, UnpackEdge(forward_hierarchy_tree[i].at(cur_node), cur_node, unpack_cache)))
                    {
                        is_first_node = false;
                        break;
                    }
                    cur_node=forward_hierarchy_tree[i].at(cur_node);
                }
                if(!is_first_node)
                {
                    blocked_cross_nodes.insert(cross_nodes_table[i * number_of_targets + idxs[j]]);
//                    SimpleLogger().Write()<<"Bloking "<<i<<" "<<idxs[j]<<". CheckIsEdgePassThrowBlocked forward.";
                    continue;
                }
                if(reached_target_phantomes.find(cur_node) != reached_target_phantomes.end())
                {
//                    SimpleLogger().Write()<<"Bloking "<<i<<" "<<idxs[j]<<". reached_target_phantomes source.";
                    continue;
                }

                cur_node = cross_nodes_table[i * number_of_targets + idxs[j]];
                while(backward_hierarchy_tree[idxs[j]].count(cur_node) && cur_node != backward_hierarchy_tree[idxs[j]].at(cur_node))
                {
                    if(SetHasIntersection(all_blocked, UnpackEdge(cur_node, backward_hierarchy_tree[idxs[j]].at(cur_node), unpack_cache)))
                    {
                        is_first_node = false;
                        break;
                    }
                    cur_node = backward_hierarchy_tree[idxs[j]].at(cur_node);
                }
                if(!is_first_node) {
//                    SimpleLogger().Write()<<"Bloking "<<i<<" "<<idxs[j]<<". CheckIsEdgePassThrowBlocked backward.";
                    continue;
                }
                if(reached_target_phantomes.find(cur_node) != reached_target_phantomes.end())
                {
//                    SimpleLogger().Write()<<"Bloking "<<i<<" "<<idxs[j]<<". reached_target_phantomes target.";
                    continue;
                }

                result_table.second[i].push_back(idxs[j]);
                if (result_table.first.first[i * number_of_targets + idxs[j]] > 0)
                    reached_target_phantomes.insert(cur_node);
//                SimpleLogger().Write()<<"Nearest for "<<i<<" is "<<idxs[j]<<" found on "<<j<<" iteration.";
            }
        }
    }

    std::set<NodeID> &UnpackEdge(NodeID from, NodeID to, std::map<std::pair<NodeID, NodeID>, std::set<NodeID>> &cache) const
    {
        auto cache_index = std::make_pair(from, to);
        auto cache_iter = cache.find(cache_index);
        if(cache_iter != cache.end())
            return cache_iter->second;
        EdgeID edge = super::facade->FindEdgeInEitherDirection(from, to);
        BOOST_ASSERT(edge != SPECIAL_EDGEID);
        EdgeData ed = super::facade->GetEdgeData(edge);
        if(ed.shortcut){
            const std::set<NodeID> &tempset1 = UnpackEdge( from, ed.id, cache);
            const std::set<NodeID> &tempset2 = UnpackEdge(ed.id,    to, cache);
            std::set_union(tempset1.begin(), tempset1.end(),
                           tempset2.begin(), tempset2.end(),
                           std::inserter(cache[cache_index], cache[cache_index].begin()));
        }
        else
            cache[cache_index].insert(from);
        return cache[cache_index];
    }

    bool SetHasIntersection(const std::set<NodeID> &set1, const std::set<NodeID> &set2) const
    {
        auto first1 = set1.begin(),
                last1 = set1.end(),
                first2 = set2.begin(),
                last2 = set2.end();
        while (first1 != last1 && first2 != last2)
        {
            if (*first1 < *first2) ++first1;
            else if (*first2 < *first1) ++first2;
            else return true;
        }
        return false;
    }

};
}
}
}

#endif
