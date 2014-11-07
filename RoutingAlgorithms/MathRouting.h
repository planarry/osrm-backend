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
#include "../DataStructures/JSONContainer.h"
#include "../DataStructures/SearchEngineData.h"
#include "../typedefs.h"

#include <boost/assert.hpp>

#include <limits>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
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
    struct SPHeapData
    {
        EdgeWeight time;
        EdgeWeight weight;
        std::vector<unsigned> settled;
        SPHeapData(const EdgeWeight time, 
                   const EdgeWeight weight, 
                   const std::vector<unsigned> &settled)
            : time(time), weight(weight), settled(settled)
        {
        }
    };
    typedef BasicRoutingInterface<DataFacadeT> super;
    typedef SearchEngineData::QueryHeap QueryHeap;
    typedef std::unordered_map<NodeID, EdgeWeight> LengthMap;
    typedef std::unordered_map<NodeID, std::vector<NodeBucket>> SearchSpaceWithBuckets;
    typedef std::vector<std::pair<NodeID, EdgeWeight>>  CrossNodesTable;
    typedef typename DataFacadeT::EdgeData EdgeData;
    SearchEngineData &engine_working_data;
    
    static const unsigned NEAREST_RADIUS = 2;

  public:
    MathRouting(DataFacadeT *facade, SearchEngineData &engine_working_data)
        : super(facade), engine_working_data(engine_working_data)
    {
    }

    ~MathRouting() {}

    std::shared_ptr<std::vector<unsigned>> operator()(const PhantomNodeArray &phantom_nodes_array, const TransportRestriction &tr, std::vector<char> &output)
        const
    {
        const unsigned number_of_locations=phantom_nodes_array.size();
        std::shared_ptr<std::vector<unsigned>> result_table =
            std::make_shared<std::vector<unsigned>>(number_of_locations,
                    std::numeric_limits<unsigned>::max());
        std::vector<EdgeWeight> time_matrix(number_of_locations*number_of_locations,
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
                                   time_matrix,
                                   tr);
            }

            ++source_id;
        }
        BOOST_ASSERT(source_id == target_id);
        
        std::vector<int> points_use_count(number_of_locations);
        
        std::vector<std::set<unsigned>> nearest_graph(number_of_locations);
        
        //i=start>|               |stoppers>>>>>>>>|  |chains>>>>>>>>>>>>>|
        std::vector<std::multimap<std::set<unsigned>, std::vector<unsigned>>> chains_with_stopers_by_start_point(number_of_locations);
        
        //i=end>>>|        one-before-end|  |start>|
        std::vector<std::multimap<unsigned, unsigned>> tails_list_ends_in_point(number_of_locations);
        
        // Building Nearest Graph
        for(unsigned i=0; i<number_of_locations; ++i)
        {
            std::vector<unsigned> sortvector(time_matrix.begin() + i * number_of_locations, 
                                             time_matrix.begin() + (i + 1) * number_of_locations);
            std::sort(sortvector.begin(), sortvector.end());
            for(unsigned j=0; j<number_of_locations; ++j)
                if(i != j && time_matrix[i * number_of_locations + j] <= sortvector[NEAREST_RADIUS])
                {
                    nearest_graph[i].insert(j);
                    SimpleLogger().Write()<<"nearest for "<<i<<" is "<<j;
                }
        }
        
        //Try Look For Chain From Every Point
        const std::set<unsigned> emptyset;
        const std::vector<std::pair<unsigned, std::multimap<std::set<unsigned>, std::vector<unsigned>>::iterator>> emptyvector;
        for(unsigned i=0; i<number_of_locations; ++i)
        {
            RecursiveLookForChain(i,
                                  emptyset, 
                                  emptyvector,
                                  chains_with_stopers_by_start_point, 
                                  tails_list_ends_in_point, 
                                  points_use_count, 
                                  nearest_graph, 
                                  time_matrix);
        }
        OutputChainsByStart(chains_with_stopers_by_start_point, output);
        return result_table;
    }
    
    void RecursiveLookForChain(const unsigned cur_point,
                               std::set<unsigned> seteled_points, 
                               std::vector<std::pair<unsigned, std::multimap<std::set<unsigned>, std::vector<unsigned>>::iterator>> chains_with_stopers_for_grow,
                               std::vector<std::multimap<std::set<unsigned>, std::vector<unsigned>>> &chains_with_stopers_by_start_point,
                               std::vector<std::multimap<unsigned, unsigned>> &tails_list_ends_in_point,
                               std::vector<int> points_use_count,
                               const std::vector<std::set<unsigned>> &nearest_graph,
                               const std::vector<EdgeWeight> &time_matrix) const
    {
        SimpleLogger().Write()<<"run from "<<cur_point<<" whith "<<seteled_points.size()<<" seteled and "<<chains_with_stopers_for_grow.size()<<" traked";
        
        std::vector<std::vector<unsigned>*> chain_from_cache;
        for(auto &chain_candidate : chains_with_stopers_by_start_point[cur_point])
        {
            
            auto cur_stopper = chain_candidate.first.begin();
            auto cur_seteled = seteled_points.begin();
            
            for(; cur_seteled != seteled_points.end(); ++cur_seteled)
            {
                if(cur_stopper != chain_candidate.first.end())
                {
                    if (*cur_stopper<*cur_seteled) break;
                    else if(*cur_seteled==*cur_stopper) ++cur_stopper;
                }
                if(std::find(chain_candidate.second.begin(),
                                  chain_candidate.second.end(),
                                  *cur_seteled) != chain_candidate.second.end())
                    break;
            }
            if(cur_seteled == seteled_points.end() && cur_stopper == chain_candidate.first.end())
                chain_from_cache.push_back(&chain_candidate.second);
        }
        
        SimpleLogger().Write()<<"chain_from_cache is "<<chain_from_cache.size();
        
        if(chain_from_cache.size() == 1)
        {
            for(auto &chain_for_grow : chains_with_stopers_for_grow)
            {
                chain_for_grow.second->second.push_back(cur_point);
                chain_for_grow.second->second.insert(chain_for_grow.second->second.end(),
                                                     chain_from_cache[0]->begin(),
                                                     chain_from_cache[0]->end());
            }
            return; //without recursion
        }
        else if(chain_from_cache.size() > 1)
        {
            for(auto &chain_for_grow : chains_with_stopers_for_grow)
                chain_for_grow.second->second.push_back(cur_point);
            return; //without recursion
        }
        
        
        
        std::set<unsigned> allowed_from_here;
        std::set<unsigned> stopped_from_here;
        auto cur_nearest = nearest_graph[cur_point].begin();
        auto cur_seteled = seteled_points.begin();
        while (cur_seteled != seteled_points.end() && cur_nearest != nearest_graph[cur_point].end())
        {
            if (*cur_nearest<*cur_seteled) 
            {
                allowed_from_here.insert(*cur_nearest);
                ++cur_nearest;
            }
            else if (*cur_seteled<*cur_nearest) ++cur_seteled;
            else 
            {
                stopped_from_here.insert(*cur_nearest);
                ++cur_nearest; 
                ++cur_seteled;
            }
        }
        if(cur_nearest != nearest_graph[cur_point].end())
            allowed_from_here.insert(cur_nearest, nearest_graph[cur_point].end());
        SimpleLogger().Write()<<"nearest is "<<nearest_graph[cur_point].size();
        SimpleLogger().Write()<<"allowed_from_here is "<<allowed_from_here.size();
        SimpleLogger().Write()<<"stopped_from_here is "<<stopped_from_here.size();
        
        
        for(auto &chain_for_grow : chains_with_stopers_for_grow) //look throw chains for grow
        {
            for(unsigned stopper : stopped_from_here) // for each stopper
                if(chain_for_grow.first != stopper //check if chain dosn't contains stopper
                    && std::find(chain_for_grow.second->second.begin(),
                                chain_for_grow.second->second.end(), 
                                stopper) == chain_for_grow.second->second.end())
                { // then update chain's stoppers
                    auto new_stoppers=chain_for_grow.second->first; //copy chain's stoppers
                    new_stoppers.insert(stopper); //add cur stopper
                    chains_with_stopers_by_start_point[chain_for_grow.first].erase(chain_for_grow.second); //delete chain from map
                    chain_for_grow.second = chains_with_stopers_by_start_point[chain_for_grow.first]
                        .emplace(new_stoppers, chain_for_grow.second->second); //insert chain with new stoppers
                }
            chain_for_grow.second->second.push_back(cur_point); //grow
        }
        
        seteled_points.insert(cur_point); // set current point as seteled
        
        if(allowed_from_here.size() == 1)
        {
            // create new chain from current point
            auto it = chains_with_stopers_by_start_point[cur_point].emplace(stopped_from_here, std::vector<unsigned>());
            //it->second.push_back(cur_point);
            chains_with_stopers_for_grow.emplace_back(cur_point, it); //track chain from current point
            RecursiveLookForChain(*allowed_from_here.begin(),
                                  seteled_points, 
                                  chains_with_stopers_for_grow,
                                  chains_with_stopers_by_start_point, 
                                  tails_list_ends_in_point, 
                                  points_use_count, 
                                  nearest_graph, 
                                  time_matrix);
        }
        else if(allowed_from_here.size() > 1)
        {
            for(unsigned next : allowed_from_here)
            {
                chains_with_stopers_for_grow.clear(); //stop traking
                // create new chain from current point
                auto it = chains_with_stopers_by_start_point[cur_point].emplace(stopped_from_here, std::vector<unsigned>());
                //it->second.push_back(cur_point);
                chains_with_stopers_for_grow.emplace_back(cur_point, it); //track chain from current point
                RecursiveLookForChain(next,
                                      seteled_points, 
                                      chains_with_stopers_for_grow,
                                      chains_with_stopers_by_start_point, 
                                      tails_list_ends_in_point, 
                                      points_use_count, 
                                      nearest_graph, 
                                      time_matrix);
            }
        }
        // else size == 0 then return without recursion
    }
    void OutputChainsByStart (const std::vector<std::multimap<std::set<unsigned>, std::vector<unsigned>>> &chains_with_stopers_by_start_point,
                              std::vector<char> &output) const
    {
        JSON::Object json_root;
        //JSON::Array chain_groups_array;
        JSON::Array chains_array;
        for(unsigned start=0; start<chains_with_stopers_by_start_point.size(); ++start)
        {
            for(const auto &chain : chains_with_stopers_by_start_point[start])
            {
                JSON::Array chain_points_array;
                chain_points_array.values.push_back(start);
                for(const unsigned point : chain.second)
                    chain_points_array.values.push_back(point);
                chains_array.values.push_back(chain_points_array);
            }
            //chain_groups_array.values.push_back(chains_array);
        }
        json_root.values["chains"] = chains_array;//chain_groups_array;
        json_root.values["n"] = chains_with_stopers_by_start_point.size();
        JSON::render(output, json_root);
    }
    
    void BuildShortestPathGaph(std::unordered_map<unsigned,std::unordered_set<long>> &start_points_for_graph,
                               std::vector<std::unordered_map<unsigned,std::unordered_set<long>>> &shortest_graph,
                               const SearchSpaceWithBuckets &backward_search_space_with_buckets,
                               const SearchSpaceWithBuckets &forward_search_space_with_buckets,
                               const std::vector<std::pair<NodeID, EdgeWeight>> &cross_nodes_table,
                               const unsigned number_of_locations) const
    {
        for(int source_id=0;source_id<number_of_locations;++source_id)
            for(int target_id=0;target_id<number_of_locations;++target_id)
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
    }

    void ForwardRoutingStep(const unsigned source_id,
                            const unsigned number_of_locations,
                            QueryHeap &query_heap,
                            const SearchSpaceWithBuckets &backward_search_space_with_buckets,
                            SearchSpaceWithBuckets &forward_search_space_with_buckets,
                            CrossNodesTable &cross_nodes_table,
                            std::vector<EdgeWeight> &time_matrix,
                            const TransportRestriction &tr) const
    {
        const NodeID node = query_heap.DeleteMin();
        const int source_distance = query_heap.GetKey(node);
        //const int source_length = length_map[node];

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
                const EdgeWeight current_distance = cross_nodes_table[index].second;
                if(!cross_nodes_table[index].first 
                    || (new_distance >= 0 && new_distance < current_distance))
                {
                    cross_nodes_table[index]=std::make_pair(node,new_distance);
                    time_matrix[index]=new_distance;
                }
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
