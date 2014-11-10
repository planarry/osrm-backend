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
#include <list>

template <class DataFacadeT> class MathRouting : public BasicRoutingInterface<DataFacadeT>
{
    typedef unsigned PointID;
    struct NodeBucket
    {
        PointID point_id; // essentially a row in the distance matrix
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
    struct WeightMatrix
    {
        const unsigned number_of_locations;
        std::vector<EdgeWeight> data;
        WeightMatrix(unsigned number_of_locations)
            : number_of_locations(number_of_locations), 
              data(number_of_locations * number_of_locations, std::numeric_limits<EdgeWeight>::max())
        {
        }
        
        EdgeWeight &at(PointID x, PointID y)
        { return data[x * number_of_locations + y]; }
        
        const EdgeWeight &at(PointID x, PointID y) const
        { return data[x * number_of_locations + y]; }
        
        std::vector<EdgeWeight>::iterator row_iter(PointID row)
        { return data.begin() + row * number_of_locations; }
        
        std::vector<EdgeWeight>::const_iterator row_iter(PointID row) const
        { return data.begin() + row * number_of_locations; }
        
    };
    struct Chain
    {
        PointID start;
        std::set<PointID> required_stoppers;
        std::vector<PointID> points_chain;
        std::set<PointID> points_set;
        bool longest;
        short level;
        EdgeWeight forward_distance;
        EdgeWeight reverse_distance;
        explicit Chain(PointID start,
                       std::set<PointID> required_stoppers,
                       bool longest, 
                       short level = 1) 
            : start(start), 
              required_stoppers(required_stoppers), 
              longest(longest), 
              level(level), 
              forward_distance(0), 
              reverse_distance(0)
        {
        }
        const PointID GetLast() const
        {
            if(points_chain.empty())
                return start;
            else return points_chain.back();
        }
        void Grow(PointID point, const WeightMatrix &weight_matrix)
        {
            forward_distance += weight_matrix.at(GetLast(), point);
            reverse_distance += weight_matrix.at(point, GetLast());
            points_chain.push_back(point);
            points_set.insert(point);
        }
        bool IsCompliant(std::set<PointID> already_attended) const
        {
            auto stopper_iter = required_stoppers.begin();
            auto attended_iter = already_attended.begin();
            
            for(; attended_iter != already_attended.end(); ++attended_iter)
            {
                if(stopper_iter != required_stoppers.end())
                {
                    if (*stopper_iter<*attended_iter) break;
                    else if(*attended_iter==*stopper_iter) ++stopper_iter;
                }
                if(points_set.find(*attended_iter) != points_set.end())
                    break;
            }
            return attended_iter == already_attended.end() && stopper_iter == required_stoppers.end();
        }
        bool operator<(const Chain &right) const
        {
            const unsigned n=std::min(points_chain.size(), right.points_chain.size());
            for(unsigned i=0; i < n; ++i)
                if(points_chain[i]<right.points_chain[i]) return true;
                else if(points_chain[i]>right.points_chain[i]) return false;
            if(points_chain.size()<right.points_chain.size()) return true;
            else return false;
        }
        bool HasPoint(PointID point) const
        {
            return start==point || points_set.find(point) != points_set.end();
        }
    };
    typedef BasicRoutingInterface<DataFacadeT> super;
    typedef SearchEngineData::QueryHeap QueryHeap;
    typedef std::unordered_map<NodeID, EdgeWeight> LengthMap;
    typedef std::unordered_map<NodeID, std::vector<NodeBucket>> SearchSpaceWithBuckets;
    typedef std::vector<std::pair<NodeID, EdgeWeight>>  CrossNodesTable;
    typedef typename DataFacadeT::EdgeData EdgeData;
    typedef typename std::vector<Chain>::iterator ChainRef;
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
        WeightMatrix time_matrix(number_of_locations);

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
        
        std::vector<std::set<PointID>> nearest_graph(number_of_locations);
        
        std::vector<std::vector<Chain>> chains_pull(number_of_locations);
        
        //i=end>>>|       one-before-end|  |start|
        std::vector<std::multimap<PointID, PointID>> tails_list(number_of_locations);
        
        // Building Nearest Graph
        for(PointID i=0; i<number_of_locations; ++i)
        {
            std::vector<EdgeWeight> sortvector(time_matrix.row_iter(i), time_matrix.row_iter(i + 1));
            std::sort(sortvector.begin(), sortvector.end());
            EdgeWeight threshold = sortvector[std::min(NEAREST_RADIUS, number_of_locations - 1)];
            for(PointID j=0; j<number_of_locations; ++j)//{SimpleLogger().Write()<<"dist from "<<i<<" to "<<j<<" is "<<time_matrix.at(i, j);
                if(i != j && time_matrix.at(i, j) <= threshold)
                {
                    nearest_graph[i].insert(j);
                    SimpleLogger().Write()<<"nearest for "<<i<<" is "<<j;
                }//}
        }
        
        //Try Look For Chain From Every Point
        const std::set<PointID> empty_attendance_set;
        const std::list<ChainRef> empty_track_list;
        for(PointID i=0; i<number_of_locations; ++i)
        {
            RecursiveLookForChain(i,
                                  empty_attendance_set, 
                                  empty_track_list,
                                  chains_pull, 
                                  tails_list, 
                                  points_use_count, 
                                  nearest_graph, 
                                  time_matrix,
                                  number_of_locations);
        }
        OutputChainsByStart(output, chains_pull);
        return result_table;
    }
    
    void RecursiveLookForChain(const PointID cur_point,
                               std::set<PointID> attendance_set, 
                               std::list<ChainRef> track_list,
                               std::vector<std::vector<Chain>> &chains_pull_ref,
                               std::vector<std::multimap<PointID, PointID>> &tails_list_ref,
                               std::vector<int> &points_use_count_ref,
                               const std::vector<std::set<PointID>> &nearest_graph,
                               const WeightMatrix &time_matrix,
                               const unsigned number_of_locations) const
    {
        SimpleLogger().Write()<<"run from "<<cur_point<<" whith "<<attendance_set.size()<<" attend and "<<track_list.size()<<" traked";
        
        std::vector<ChainRef> chain_from_cache;
        for(ChainRef chain_candidate = chains_pull_ref[cur_point].begin(); 
            chain_candidate < chains_pull_ref[cur_point].end();
            ++chain_candidate)
            if(chain_candidate->IsCompliant(attendance_set))
                chain_from_cache.push_back(chain_candidate);
        
        SimpleLogger().Write()<<"chain_from_cache is "<<chain_from_cache.size();
        
        if(chain_from_cache.size() == 1)
        {
            chain_from_cache.front()->longest = false;
            for(ChainRef chain_for_grow : track_list)
            {
                chain_for_grow->Grow(cur_point, time_matrix);
                for(PointID point : chain_from_cache.front()->points_chain)
                    chain_for_grow->Grow(point, time_matrix);
            }
            return; //without recursion
        }
        else if(chain_from_cache.size() > 1)
        {
            for(ChainRef chain_for_grow : track_list)
                chain_for_grow->Grow(cur_point, time_matrix);
            return; //without recursion
        }
        
        
        
        std::set<PointID> allowed_from_here;
        std::set<PointID> stopped_from_here;
        auto nearest_iter = nearest_graph[cur_point].begin();
        auto attended_iter = attendance_set.begin();
        while (attended_iter != attendance_set.end() && nearest_iter != nearest_graph[cur_point].end())
        {
            if (*nearest_iter<*attended_iter) 
            {
                allowed_from_here.insert(*nearest_iter);
                ++nearest_iter;
            }
            else if (*attended_iter<*nearest_iter) ++attended_iter;
            else 
            {
                stopped_from_here.insert(*nearest_iter);
                ++nearest_iter; 
                ++attended_iter;
            }
        }
        if(nearest_iter != nearest_graph[cur_point].end())
            allowed_from_here.insert(nearest_iter, nearest_graph[cur_point].end());
        SimpleLogger().Write()<<"nearest is "<<nearest_graph[cur_point].size();
        SimpleLogger().Write()<<"allowed_from_here is "<<allowed_from_here.size();
        SimpleLogger().Write()<<"stopped_from_here is "<<stopped_from_here.size();
        
        //look throw traked chains for grow
        for(ChainRef chain_for_grow : track_list) 
        {
            //update stoppers
            for(PointID stopper : stopped_from_here)
                if(!chain_for_grow->HasPoint(stopper))
                    chain_for_grow->required_stoppers.insert(stopper);
            SimpleLogger().Write()<<"grow";
            chain_for_grow->Grow(cur_point, time_matrix);
        }
        
        attendance_set.insert(cur_point); // set current point as attended
        
        if(allowed_from_here.size() == 1)
        {
            // create new chain from current point
            chains_pull_ref[cur_point].emplace_back(cur_point, stopped_from_here, false);
            track_list.emplace_back(chains_pull_ref[cur_point].end() - 1); //track chain from current point
            RecursiveLookForChain(*allowed_from_here.begin(),
                                  attendance_set, 
                                  track_list,
                                  chains_pull_ref, 
                                  tails_list_ref, 
                                  points_use_count_ref, 
                                  nearest_graph, 
                                  time_matrix,
                                  number_of_locations);
        }
        else if(allowed_from_here.size() > 1)
        {
            for(PointID next : allowed_from_here)
            {
                track_list.clear(); //stop traking
                // create new chain from current point
                chains_pull_ref[cur_point].emplace_back(cur_point, stopped_from_here, true);
                track_list.emplace_back(chains_pull_ref[cur_point].end() - 1); //track chain from current point
                RecursiveLookForChain(next,
                                      attendance_set, 
                                      track_list,
                                      chains_pull_ref, 
                                      tails_list_ref, 
                                      points_use_count_ref, 
                                      nearest_graph, 
                                      time_matrix,
                                      number_of_locations);
            }
        }
        // else size == 0 then return without recursion
    }
    void OutputChainsByStart (std::vector<char> &output,
                              const std::vector<std::vector<Chain>> &chains_pull) const
    {
        JSON::Object json_root;
        //JSON::Array chain_groups_array;
        JSON::Array chains_array;
        for(PointID start=0; start<chains_pull.size(); ++start)
        {
            for(const Chain &chain : chains_pull[start])
                if(chain.longest)
                {
                    JSON::Array chain_points_array;
                    chain_points_array.values.push_back(start);
                    for(const PointID point : chain.points_chain)
                        chain_points_array.values.push_back(point);
                    chains_array.values.push_back(chain_points_array);
                }
            //chain_groups_array.values.push_back(chains_array);
        }
        json_root.values["chains"] = chains_array;//chain_groups_array;
        json_root.values["n"] = chains_pull.size();
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
                            WeightMatrix &time_matrix,
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
                if(new_distance >= 0 && (!cross_nodes_table[index].first || new_distance < current_distance))
                {
                    cross_nodes_table[index] = std::make_pair(node,new_distance);
                    time_matrix.at(source_id, target_id) = new_distance;
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
