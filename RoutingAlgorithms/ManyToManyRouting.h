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

#ifndef MANY_TO_MANY_ROUTING_H
#define MANY_TO_MANY_ROUTING_H

#include "BasicRoutingInterface.h"
#include "../DataStructures/SearchEngineData.h"
#include "../typedefs.h"

#include <boost/assert.hpp>

#include <limits>
#include <memory>
#include <unordered_map>
#include <vector>

template <class DataFacadeT> class ManyToManyRouting : public BasicRoutingInterface<DataFacadeT>
{
    typedef BasicRoutingInterface<DataFacadeT> super;
    typedef SearchEngineData::QueryHeap QueryHeap;
    SearchEngineData &engine_working_data;
    std::shared_ptr<std::unordered_map<NodeID, EdgeWeight>> length_map;

    struct NodeBucket
    {
        unsigned target_id; // essentially a row in the distance matrix
        EdgeWeight distance;
        EdgeWeight length;
        NodeBucket(const unsigned target_id, const EdgeWeight distance, const EdgeWeight length)
            : target_id(target_id), distance(distance), length(length)
        {
        }
    };
    typedef std::unordered_map<NodeID, std::vector<NodeBucket>> SearchSpaceWithBuckets;

  public:
    ManyToManyRouting(DataFacadeT *facade, SearchEngineData &engine_working_data)
        : super(facade), engine_working_data(engine_working_data),
        length_map(std::make_shared<std::unordered_map<NodeID, EdgeWeight>>())
    {
    }

    ~ManyToManyRouting() {}

    std::shared_ptr<std::vector<EdgeWeight>> operator()(const PhantomNodeArray &phantom_nodes_array)
        const
    {
        const unsigned number_of_locations = static_cast<unsigned>(phantom_nodes_array.size());
        std::shared_ptr<std::vector<EdgeWeight>> result_table =
            std::make_shared<std::vector<EdgeWeight>>(number_of_locations * number_of_locations * 2,
                                                      std::numeric_limits<EdgeWeight>::max());

        engine_working_data.InitializeOrClearFirstThreadLocalStorage(
            super::facade->GetNumberOfNodes());

        QueryHeap &query_heap = *(engine_working_data.forwardHeap);

        SearchSpaceWithBuckets search_space_with_buckets;

        unsigned target_id = 0;
        for (const std::vector<PhantomNode> &phantom_node_vector : phantom_nodes_array)
        {
            query_heap.Clear();
            length_map->clear();
            // insert target(s) at distance 0

            for (const PhantomNode &phantom_node : phantom_node_vector)
            {
                if (SPECIAL_NODEID != phantom_node.forward_node_id)
                {
                    query_heap.Insert(phantom_node.forward_node_id,
                                      phantom_node.GetForwardWeightPlusOffset(),
                                      phantom_node.forward_node_id);
                    (*length_map)[phantom_node.forward_node_id]=phantom_node.GetForwardLength();
                }
                if (SPECIAL_NODEID != phantom_node.reverse_node_id)
                {
                    query_heap.Insert(phantom_node.reverse_node_id,
                                      phantom_node.GetReverseWeightPlusOffset(),
                                      phantom_node.reverse_node_id);
                    (*length_map)[phantom_node.reverse_node_id]=phantom_node.GetReverseLength();
                }
            }

            // explore search space
            while (!query_heap.Empty())
            {
                BackwardRoutingStep(target_id, query_heap, search_space_with_buckets);
            }
            ++target_id;
        }

        // for each source do forward search
        unsigned source_id = 0;
        for (const std::vector<PhantomNode> &phantom_node_vector : phantom_nodes_array)
        {
            query_heap.Clear();
            length_map->clear();
            for (const PhantomNode &phantom_node : phantom_node_vector)
            {
                // insert sources at distance 0
                if (SPECIAL_NODEID != phantom_node.forward_node_id)
                {
                    query_heap.Insert(phantom_node.forward_node_id,
                                      -phantom_node.GetForwardWeightPlusOffset(),
                                      phantom_node.forward_node_id);
                    (*length_map)[phantom_node.forward_node_id]=-phantom_node.GetForwardLength();
                }
                if (SPECIAL_NODEID != phantom_node.reverse_node_id)
                {
                    query_heap.Insert(phantom_node.reverse_node_id,
                                      -phantom_node.GetReverseWeightPlusOffset(),
                                      phantom_node.reverse_node_id);
                    (*length_map)[phantom_node.reverse_node_id]=-phantom_node.GetReverseLength();
                }
            }

            // explore search space
            while (!query_heap.Empty())
            {
                ForwardRoutingStep(source_id,
                                   number_of_locations,
                                   query_heap,
                                   search_space_with_buckets,
                                   result_table);
            }

            ++source_id;
        }
        BOOST_ASSERT(source_id == target_id);
        return result_table;
    }

    void ForwardRoutingStep(const unsigned source_id,
                            const unsigned number_of_locations,
                            QueryHeap &query_heap,
                            const SearchSpaceWithBuckets &search_space_with_buckets,
                            std::shared_ptr<std::vector<EdgeWeight>> result_table) const
    {
        const NodeID node = query_heap.DeleteMin();
        const int source_distance = query_heap.GetKey(node);
        const int source_length = (*length_map)[node];

        // check if each encountered node has an entry
        const auto bucket_iterator = search_space_with_buckets.find(node);
        // iterate bucket if there exists one
        if (bucket_iterator != search_space_with_buckets.end())
        {
            const std::vector<NodeBucket> &bucket_list = bucket_iterator->second;
            for (const NodeBucket &current_bucket : bucket_list)
            {
                // get target id from bucket entry
                const unsigned target_id = current_bucket.target_id;
                const int target_distance = current_bucket.distance;
                const EdgeWeight current_distance =
                    (*result_table)[source_id * number_of_locations + target_id];
                // check if new distance is better
                const EdgeWeight new_distance = source_distance + target_distance;
                if (new_distance >= 0 && new_distance < current_distance)
                {
                    (*result_table)[source_id * number_of_locations + target_id] =
                        (source_distance + target_distance);
                    (*result_table)[source_id * number_of_locations + target_id + number_of_locations*number_of_locations] =
                        (source_length + current_bucket.length);
                }
            }
        }
        if (StallAtNode<true>(node, source_distance, query_heap))
        {
            return;
        }
        RelaxOutgoingEdges<true>(node, source_distance, query_heap);
    }

    void BackwardRoutingStep(const unsigned target_id,
                             QueryHeap &query_heap,
                             SearchSpaceWithBuckets &search_space_with_buckets) const
    {
        const NodeID node = query_heap.DeleteMin();
        const int target_distance = query_heap.GetKey(node);

        // store settled nodes in search space bucket
        search_space_with_buckets[node].emplace_back(target_id, target_distance, (*length_map)[node]);

        if (StallAtNode<false>(node, target_distance, query_heap))
        {
            return;
        }

        RelaxOutgoingEdges<false>(node, target_distance, query_heap);
    }

    template <bool forward_direction>
    inline void
    RelaxOutgoingEdges(const NodeID node, const EdgeWeight distance, QueryHeap &query_heap) const
    {
        for (auto edge : super::facade->GetAdjacentEdgeRange(node))
        {
            const auto &data = super::facade->GetEdgeData(edge);
            const bool direction_flag = (forward_direction ? data.forward : data.backward);
            if (direction_flag)
            {
                const NodeID to = super::facade->GetTarget(edge);
                const int edge_weight = data.distance;

                BOOST_ASSERT_MSG(edge_weight > 0, "edge_weight invalid");
                const int to_distance = distance + edge_weight;
                const int to_length = (*length_map)[node] + data.length;

                // New Node discovered -> Add to Heap + Node Info Storage
                if (!query_heap.WasInserted(to))
                {
                    query_heap.Insert(to, to_distance, node);
                    (*length_map)[to]=to_length;
                }
                // Found a shorter Path -> Update distance
                else if (to_distance < query_heap.GetKey(to))
                {
                    // new parent
                    query_heap.GetData(to).parent = node;
                    query_heap.DecreaseKey(to, to_distance);
                    (*length_map)[to]=to_length;
                }
            }
        }
    }

    // Stalling
    template <bool forward_direction>
    inline bool StallAtNode(const NodeID node, const EdgeWeight distance, QueryHeap &query_heap)
        const
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
