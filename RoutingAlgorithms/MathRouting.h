/*

Copyright (c) 2014, Project OSRM, Dennis Luxen, others
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list
of conditions and the following disclaimer.
Redistributions in binary form must reproduce the above copyright notice, this
list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef MATH_ROUTING_H
#define MATH_ROUTING_H

#include "BasicRoutingInterface.h"
#include "../DataStructures/SearchEngineData.h"
#include "../typedefs.h"

#include <boost/assert.hpp>

#include <limits>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

template <class DataFacadeT> class MathRouting : public BasicRoutingInterface<DataFacadeT>
{
    struct NodeBucket
    {
        unsigned point_id; // essentially a row in the distance matrix
        EdgeWeight distance;
        NodeID parent;
        NodeBucket(const unsigned point_id, const EdgeWeight distance, const NodeID parent)
            : point_id(point_id), distance(distance), parent(parent)
        {
        }
    };
    typedef BasicRoutingInterface<DataFacadeT> super;
    typedef SearchEngineData::QueryHeap QueryHeap;
    typedef std::unordered_map<NodeID, EdgeWeight> LengthMap;
    typedef std::unordered_map<NodeID, std::vector<NodeBucket>> SearchSpaceWithBuckets;
    typedef std::vector<std::pair<NodeID, EdgeWeight>>  CrossNodesTable;
    SearchEngineData &engine_working_data;

  public:
    MathRouting(DataFacadeT *facade, SearchEngineData &engine_working_data)
        : super(facade), engine_working_data(engine_working_data)
    {
    }

    ~MathRouting() {}

    std::shared_ptr<std::vector<unsigned>> operator()(const PhantomNodeArray &phantom_nodes_array, const TransportRestriction &tr)
        const
    {
        const unsigned number_of_locations=phantom_nodes_array.size();
        std::shared_ptr<std::vector<unsigned>> result_table =
            std::make_shared<std::vector<unsigned>>(number_of_locations,
                    std::numeric_limits<unsigned>::max());

        engine_working_data.InitializeOrClearFirstThreadLocalStorage(
            super::facade->GetNumberOfNodes());

        QueryHeap &query_heap = *(engine_working_data.forwardHeap);

        SearchSpaceWithBuckets backward_search_space_with_buckets, forward_search_space_with_buckets;
        std::vector<std::pair<NodeID, EdgeWeight>> cross_nodes_table(number_of_locations*number_of_locations);

        unsigned target_id = 0;
        for (const std::vector<PhantomNode> &phantom_node_vector : phantom_nodes_array)
        {
            query_heap.Clear();
            // insert target(s) at distance 0

            for (const PhantomNode &phantom_node : phantom_node_vector)
            {
                if (SPECIAL_NODEID != phantom_node.forward_node_id)
                {
                    query_heap.Insert(phantom_node.forward_node_id,
                                      phantom_node.GetForwardWeightPlusOffset(),
                                      phantom_node.forward_node_id);
                }
                if (SPECIAL_NODEID != phantom_node.reverse_node_id)
                {
                    query_heap.Insert(phantom_node.reverse_node_id,
                                      phantom_node.GetReverseWeightPlusOffset(),
                                      phantom_node.reverse_node_id);
                }
            }

            // explore search space
            while (!query_heap.Empty())
            {
                BackwardRoutingStep(target_id, query_heap, backward_search_space_with_buckets, tr);
            }
            ++target_id;
        }

        // for each source do forward search
        unsigned source_id = 0;
        for (const std::vector<PhantomNode> &phantom_node_vector : phantom_nodes_array)
        {
            query_heap.Clear();
            for (const PhantomNode &phantom_node : phantom_node_vector)
            {
                // insert sources at distance 0
                if (SPECIAL_NODEID != phantom_node.forward_node_id)
                {
                    query_heap.Insert(phantom_node.forward_node_id,
                                      -phantom_node.GetForwardWeightPlusOffset(),
                                      phantom_node.forward_node_id);
                }
                if (SPECIAL_NODEID != phantom_node.reverse_node_id)
                {
                    query_heap.Insert(phantom_node.reverse_node_id,
                                      -phantom_node.GetReverseWeightPlusOffset(),
                                      phantom_node.reverse_node_id);
                }
            }

            // explore search space
            while (!query_heap.Empty())
            {
                ForwardRoutingStep(source_id,
                                   number_of_locations,
                                   query_heap, 
                                   backward_search_space_with_buckets,
                                   forward_search_space_with_buckets,
                                   cross_nodes_table,
                                   tr);
            }

            ++source_id;
        }
        BOOST_ASSERT(source_id == target_id);
        
        std::unordered_map<unsigned,std::unordered_set<long>> start_points_for_graph;
        std::vector<std::unordered_map<unsigned,std::unordered_set<long>>> shortest_graph(number_of_locations);
        for(source_id=0;source_id<number_of_locations;++source_id)
            for(target_id=0;target_id<number_of_locations;++target_id)
            {
                NodeID cur_node=cross_nodes_table[source_id * number_of_locations + target_id].first;
                bool process=true;
                while(process)
                    for (const NodeBucket &current_bucket : backward_search_space_with_buckets[cur_node])
                        if(current_bucket.point_id==target_id)
                        {
                            if(cur_node!=current_bucket.parent)
                            {
                                shortest_graph[source_id][cur_node].insert(current_bucket.parent);
                                cur_node==current_bucket.parent;
                            }
                            else 
                            {
                                shortest_graph[source_id][cur_node].insert(-target_id);
                                process=false;
                            }
                            break;
                        }
                cur_node=cross_nodes_table[source_id * number_of_locations + target_id].first;
                process=true;
                while(process)
                    for (const NodeBucket &current_bucket : forward_search_space_with_buckets[cur_node])
                        if(current_bucket.point_id==source_id)
                        {
                            if(cur_node!=current_bucket.parent)
                            {
                                shortest_graph[source_id][current_bucket.parent].insert(cur_node);
                                cur_node==current_bucket.parent;
                            }
                            else 
                            {
                                process=false;
                                start_points_for_graph[source_id].insert(cur_node);
                            }
                            break;
                        }
            }
        
        int current_id=0;
        
        
        
        return result_table;
    }

    void ForwardRoutingStep(const unsigned source_id,
                            const unsigned number_of_locations,
                            QueryHeap &query_heap,
                            LengthMap &length_map,
                            const SearchSpaceWithBuckets &backward_search_space_with_buckets,
                            SearchSpaceWithBuckets &forward_search_space_with_buckets,
                            CrossNodesTable &cross_nodes_table,
                            const TransportRestriction &tr) const
    {
        const NodeID node = query_heap.DeleteMin();
        const int source_distance = query_heap.GetKey(node);
        const int source_length = length_map[node];

        forward_search_space_with_buckets[node].emplace_back(source_id, source_distance, query_heap.GetData(node).parent);

        
        // check if each encountered node has an entry
        const auto bucket_iterator = backward_search_space_with_buckets.find(node);
        // iterate bucket if there exists one
        if (bucket_iterator != backward_search_space_with_buckets.end())
        {
            const std::vector<NodeBucket> &bucket_list = bucket_iterator->second;
            for (const NodeBucket &current_bucket : bucket_list)
            {
                // get target id from bucket entry
                const unsigned target_id = current_bucket.point_id;
                const int target_distance = current_bucket.distance;
                const EdgeWeight new_distance = source_distance + target_distance;
                const unsigned index = source_id * number_of_locations + target_id;
                if(cross_nodes_table[index].first) {
                    const EdgeWeight current_distance = cross_nodes_table[index].second;
                    if (new_distance >= 0 && new_distance < current_distance)
                        cross_nodes_table[index]=std::make_pair(node,new_distance);
                }
                else cross_nodes_table[index]=std::make_pair(node,new_distance);
            }
        }
        if (StallAtNode<true>(node, source_distance, query_heap, tr))
        {
            return;
        }
        RelaxOutgoingEdges<true>(node, source_distance, query_heap, tr);
    }

    void BackwardRoutingStep(const unsigned target_id,
                             QueryHeap &query_heap,
                             SearchSpaceWithBuckets &backward_search_space_with_buckets, 
                             const TransportRestriction &tr) const
    {
        const NodeID node = query_heap.DeleteMin();
        const int target_distance = query_heap.GetKey(node);

        // store settled nodes in search space bucket
        backward_search_space_with_buckets[node].emplace_back(target_id, target_distance, query_heap.GetData(node).parent);

        if (StallAtNode<false>(node, target_distance, query_heap, tr))
        {
            return;
        }

        RelaxOutgoingEdges<false>(node, target_distance, query_heap, tr);
    }

    template <bool forward_direction>
    inline void
    RelaxOutgoingEdges(const NodeID node, const EdgeWeight distance, QueryHeap &query_heap, const TransportRestriction &tr) const
    {
        for (auto edge : super::facade->GetAdjacentEdgeRange(node))
        {
            const auto &data = super::facade->GetEdgeData(edge);
            if(tr.IsEdgeRestricted(data)) continue;
            const bool direction_flag = (forward_direction ? data.forward : data.backward);
            if (direction_flag)
            {
                const NodeID to = super::facade->GetTarget(edge);
                const int edge_weight = data.distance;

                BOOST_ASSERT_MSG(edge_weight > 0, "edge_weight invalid");
                const int to_distance = distance + edge_weight;

                // New Node discovered -> Add to Heap + Node Info Storage
                if (!query_heap.WasInserted(to))
                {
                    query_heap.Insert(to, to_distance, node);
                }
                // Found a shorter Path -> Update distance
                else if (to_distance < query_heap.GetKey(to))
                {
                    // new parent
                    query_heap.GetData(to).parent = node;
                    query_heap.DecreaseKey(to, to_distance);
                }
            }
        }
    }

    // Stalling
    template <bool forward_direction>
    inline bool StallAtNode(const NodeID node, const EdgeWeight distance, QueryHeap &query_heap, const TransportRestriction &tr)
        const
    {
        for (auto edge : super::facade->GetAdjacentEdgeRange(node))
        {
            const auto &data = super::facade->GetEdgeData(edge);
            if(tr.IsEdgeRestricted(data)) continue;
            const bool reverse_flag = ((!forward_direction) ? data.forward : data.backward);
            if (reverse_flag)
            {
                const NodeID to = super::facade->GetTarget(edge);
                const int edge_weight = data.distance;
                
                BOOST_ASSERT_MSG(edge_weight > 0, "edge_weight invalid");
                if (query_heap.WasInserted(to))
                {
                    if (query_heap.GetKey(to) + edge_weight < distance)
                    {
                        return true;
                    }
                }
            }
        }
        return false;
    }
};
#endif
