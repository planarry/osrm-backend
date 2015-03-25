#ifndef MATH_ROUTING_H
#define MATH_ROUTING_H

#include "BasicRoutingInterface.h"
#include "../DataStructures/JSONContainer.h"
#include "../DataStructures/SearchEngineData.h"
#include "../logistic/GraphLogistic.h"
#include "../typedefs.h"
#include "../Util/TimingUtil.h"

#include <boost/assert.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <limits>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
#include <vector>
#include <list>

namespace ublas = boost::numeric::ublas;

template <class DataFacadeT> class MathRouting : public BasicRoutingInterface<DataFacadeT>
{
    typedef unsigned PointID;
    struct NodeBucket
    {
        PointID target_id; // essentially a row in the distance matrix
        EdgeWeight distance;
        EdgeWeight length;
        NodeID parent;
        NodeBucket(const unsigned target_id, const EdgeWeight distance, const EdgeWeight length, const NodeID parent)
            : target_id(target_id), distance(distance), length(length), parent(parent)
        {
        }
    };
    class DistanceComparator
    {
        const PointID from;
        const ublas::matrix<EdgeWeight> &matrix; 
    public:
        DistanceComparator(const PointID from, const ublas::matrix<EdgeWeight> &matrix)
            : from(from), matrix(matrix)
        {
        }
        bool operator()(PointID a, PointID b) const 
        {
            EdgeWeight a_dist, b_dist;
            if(a < 0) a_dist = matrix(1 - a, from);
            else a_dist = matrix(from, a);
            if(b < 0) b_dist = matrix(1 - b, from);
            else b_dist = matrix(from, b);
            return a_dist < b_dist;
        }
    };
    typedef BasicRoutingInterface<DataFacadeT> super;
    typedef SearchEngineData::QueryHeap QueryHeap;
    typedef std::unordered_map<NodeID, EdgeWeight> LengthMap;
    typedef std::unordered_map<NodeID, std::vector<NodeBucket>> SearchSpaceWithBuckets;
    typedef std::vector<std::unordered_map<NodeID, NodeID>> HierarchyTree;
    //typedef std::vector<NodeID>  CrossNodesTable;
    typedef typename DataFacadeT::EdgeData EdgeData;
    //typedef typename std::vector<Chain>::iterator ChainRef;
    SearchEngineData &engine_working_data;
    
    //static const unsigned NEAREST_RADIUS = 2;

  public:
    MathRouting(DataFacadeT *facade, SearchEngineData &engine_working_data)
        : super(facade), engine_working_data(engine_working_data)
    {
    }

    ~MathRouting() {}

    void operator()(const PhantomNodeArray &phantom_nodes_array, 
                    const std::vector<FixedPointCoordinate> &coordinates,  
                    const TransportRestriction &tr, 
                    std::vector<char> &output)
        const
    {
        //TIMER_START(process);
        const unsigned n = phantom_nodes_array.size();
        ublas::matrix<NodeID> cross_nodes_table(n, n);
        ublas::matrix<EdgeWeight> time_matrix(n, n),
                                  length_matrix(n, n);
        //std::fill(time_matrix.begin1(),time_matrix.end1(),std::numeric_limits<EdgeWeight>::max());
        //std::fill(length_matrix.begin1(),length_matrix.end1(),std::numeric_limits<EdgeWeight>::max());
        time_matrix = length_matrix = ublas::scalar_matrix<EdgeWeight>(n, n, std::numeric_limits<EdgeWeight>::max());
        
        engine_working_data.InitializeOrClearFirstThreadLocalStorage(
            super::facade->GetNumberOfNodes());

        QueryHeap &query_heap = *(engine_working_data.forwardHeap);
        LengthMap length_map;

        SearchSpaceWithBuckets search_space_with_buckets;
        HierarchyTree backward_hierarchy_tree(n), forward_hierarchy_tree(n);
        
        std::vector<std::set<NodeID>> phantomes_for_point(n);
        std::set<NodeID> all_phantomes;
        
        unsigned target_id = 0;
        for (const std::vector<PhantomNode> &phantom_node_vector : phantom_nodes_array)
        {
            query_heap.Clear();
            length_map.clear();
            // insert target(s) at distance 0

            for (const PhantomNode &phantom_node : phantom_node_vector)
            {
                if (SPECIAL_NODEID != phantom_node.forward_node_id)
                {
                    phantomes_for_point[target_id].insert(phantom_node.forward_node_id);
                    all_phantomes.insert(phantom_node.forward_node_id);
                    query_heap.Insert(phantom_node.forward_node_id,
                                      phantom_node.GetForwardWeightPlusOffset(),
                                      phantom_node.forward_node_id);
                    length_map[phantom_node.forward_node_id] = phantom_node.GetForwardLength();
                }
                if (SPECIAL_NODEID != phantom_node.reverse_node_id)
                {
                    phantomes_for_point[target_id].insert(phantom_node.reverse_node_id);
                    all_phantomes.insert(phantom_node.reverse_node_id);
                    query_heap.Insert(phantom_node.reverse_node_id,
                                      phantom_node.GetReverseWeightPlusOffset(),
                                      phantom_node.reverse_node_id);
                    length_map[phantom_node.reverse_node_id] = phantom_node.GetReverseLength();
                }
            }

            // explore search space
            while (!query_heap.Empty())
            {
                BackwardRoutingStep(target_id, 
                                    query_heap, 
                                    length_map, 
                                    search_space_with_buckets, 
                                    backward_hierarchy_tree, 
                                    tr);
            }
            ++target_id;
        }

        // for each source do forward search
        unsigned source_id = 0;
        for (const std::vector<PhantomNode> &phantom_node_vector : phantom_nodes_array)
        {
            query_heap.Clear();
            length_map.clear();
            for (const PhantomNode &phantom_node : phantom_node_vector)
            {
                // insert sources at distance 0
                if (SPECIAL_NODEID != phantom_node.forward_node_id)
                {
                    query_heap.Insert(phantom_node.forward_node_id,
                                      -phantom_node.GetForwardWeightPlusOffset(),
                                      phantom_node.forward_node_id);
                    length_map[phantom_node.forward_node_id] = phantom_node.GetForwardLength();
                }
                if (SPECIAL_NODEID != phantom_node.reverse_node_id)
                {
                    query_heap.Insert(phantom_node.reverse_node_id,
                                      -phantom_node.GetReverseWeightPlusOffset(),
                                      phantom_node.reverse_node_id);
                    length_map[phantom_node.reverse_node_id] = phantom_node.GetReverseLength();
                }
            }

            // explore search space
            while (!query_heap.Empty())
            {
                ForwardRoutingStep(source_id,
                                   n,
                                   query_heap, 
                                   length_map,
                                   search_space_with_buckets,
                                   forward_hierarchy_tree,
                                   backward_hierarchy_tree,
                                   time_matrix,
                                   length_matrix,
                                   cross_nodes_table,
                                   tr);
            }

            ++source_id;
        }
        BOOST_ASSERT(source_id == target_id);
        
        
        
        std::vector<std::list<PointID>> nearest_forward_graph(n), nearest_reverse_graph(n);
        //TIMER_START(solution);
        TIMER_START(graph);
        BuildFullGraph(nearest_forward_graph, 
                       nearest_reverse_graph, 
                       time_matrix,
                       cross_nodes_table, 
                       forward_hierarchy_tree, 
                       backward_hierarchy_tree,
                       phantomes_for_point,
                       all_phantomes, 
                       n);
        TIMER_STOP(graph);
        SimpleLogger().Write() << "Ful graph " << TIMER_SEC(graph) << "s";
        
        JSON::Object json_object;
        JSON::Array json_matrix_time;
        JSON::Array json_matrix_length;
        JSON::Array json_matrix_graph;
        for (unsigned row = 0; row < n; ++row)
        {
            JSON::Array json_row_time;
            JSON::Array json_row_length;
            JSON::Array json_row_graph;
            for (unsigned col = 0; col < n; ++col)
            {
                json_row_time.values.push_back(int(time_matrix(row,col)/10));
                json_row_length.values.push_back(int(length_matrix(row,col)/10));
            }
            json_row_graph.values.insert(json_row_graph.values.end(), 
                                         nearest_forward_graph[row].begin(), 
                                         nearest_forward_graph[row].end());
            json_matrix_time.values.push_back(json_row_time);
            json_matrix_length.values.push_back(json_row_length);
            json_matrix_graph.values.push_back(json_row_graph);
        }
        json_object.values["time_table"] = json_matrix_time;
        json_object.values["length_table"] = json_matrix_length;
        json_object.values["graph_table"] = json_matrix_graph;
        JSON::render(output, json_object);
        /*TIMER_START(math);
        GraphLogistic math(n, 
                           time_matrix, 
                           length_matrix, 
                           nearest_forward_graph, 
                           coordinates);
        math.run();
        TIMER_STOP(math);
        TIMER_STOP(solution);
        math.render(output, TIMER_SEC(solution));
        SimpleLogger().Write() << "Math takes " << TIMER_SEC(math) << "s";
        SimpleLogger().Write() << "Common solution takes " << TIMER_SEC(solution) << "s";*/
        
    }
    
    void BuildFullGraph(std::vector<std::list<PointID>> &nearest_forward_graph,
                        std::vector<std::list<PointID>> &nearest_reverse_graph,
                        const ublas::matrix<EdgeWeight> &time_matrix, 
                        const ublas::matrix<NodeID> &cross_nodes_table, 
                        const HierarchyTree &forward_hierarchy_tree, 
                        const HierarchyTree &backward_hierarchy_tree,
                        const std::vector<std::set<NodeID>> &phantomes_for_point,
                        const std::set<NodeID> &all_phantomes,
                        const unsigned n) const
    {
        std::map<std::pair<NodeID, NodeID>, std::set<NodeID>> unpack_cache;
        for(PointID i = 0; i < n; ++i)
        {
            std::vector<PointID> idxs;
            for(PointID j = 0; j < n; ++j)
                if(i != j) idxs.push_back(j);
            std::sort(idxs.begin(), idxs.end(), [i, &time_matrix](PointID l, PointID r){
                return time_matrix(i, l)<time_matrix(i, r);
            });
            /*TIMER_START(source);
            EdgeWeight threshold = std::numeric_limits<EdgeWeight>::max()-1;//idxs[n-2]/2;
            SimpleLogger().Write()<<"Initial threshold for "<<i<<"="<<threshold;*/
            std::set<NodeID> blocked_cross_nodes, 
                             all_blocked, 
                             reached_target_phantomes;
            for(PointID j = 0; j < std::min(n - 1, 100u); ++j)
            {
                TIMER_START(target);
                NodeID cur_node = cross_nodes_table(i, idxs[j]);
                if(blocked_cross_nodes.find(cur_node) != blocked_cross_nodes.end())
                {
                    TIMER_STOP(target);
                    SimpleLogger().Write()<<"Bloking "<<i<<" "<<idxs[j]<<". blocked_cross_nodes.\t after "<<TIMER_SEC(target)<<"s"; 
                    continue;
                }
                /*if(SetHasIntesection(reached_target_phantomes, phantomes_for_point[idxs[j]]))
                {
                    TIMER_STOP(target);
                    SimpleLogger().Write()<<"Bloking "<<i<<" "<<idxs[j]<<". target phantom_node already reached.\t after "<<TIMER_SEC(target)<<"s";
                    continue;
                }*/
                                
                all_blocked.clear();
                std::set<NodeID> tempset;        
                std::set_union(phantomes_for_point[i].begin(), phantomes_for_point[i].end(),
                               phantomes_for_point[idxs[j]].begin(), phantomes_for_point[idxs[j]].end(),
                               std::inserter(tempset, tempset.begin()));
                std::set_difference(all_phantomes.begin(), all_phantomes.end(),
                                    tempset.begin(), tempset.end(),
                                    std::inserter(all_blocked, all_blocked.begin()));
                all_blocked.insert(blocked_cross_nodes.begin(), blocked_cross_nodes.end());
                bool is_first_node = true;
                
                while(forward_hierarchy_tree[i].count(cur_node) && cur_node != forward_hierarchy_tree[i].at(cur_node))
                {
                    if(SetHasIntesection(all_blocked, UnpackEdge(forward_hierarchy_tree[i].at(cur_node), cur_node, unpack_cache)))
                    {
                        is_first_node=false;
                        break;
                    }
                    cur_node=forward_hierarchy_tree[i].at(cur_node);
                }
                if(!is_first_node)
                {
                    blocked_cross_nodes.insert(cross_nodes_table(i, idxs[j]));
                    TIMER_STOP(target);
                    SimpleLogger().Write()<<"Bloking "<<i<<" "<<idxs[j]<<". CheckIsEdgePassThrowBlocked forward.\t after "<<TIMER_SEC(target)<<"s";
                    continue;
                }
                if(reached_target_phantomes.find(cur_node) != reached_target_phantomes.end())
                {
                    TIMER_STOP(target);
                    SimpleLogger().Write()<<"Bloking "<<i<<" "<<idxs[j]<<". reached_target_phantomes source.\t after "<<TIMER_SEC(target)<<"s";
                    continue;
                }
                
                cur_node = cross_nodes_table(i, idxs[j]);
                while(backward_hierarchy_tree[idxs[j]].count(cur_node) && cur_node != backward_hierarchy_tree[idxs[j]].at(cur_node))
                {
                    if(SetHasIntesection(all_blocked, UnpackEdge(cur_node, backward_hierarchy_tree[idxs[j]].at(cur_node), unpack_cache)))
                    {
                        is_first_node=false;
                        break;
                    }
                    cur_node=backward_hierarchy_tree[idxs[j]].at(cur_node);
                }
                if(!is_first_node) {
                    TIMER_STOP(target);
                    SimpleLogger().Write()<<"Bloking "<<i<<" "<<idxs[j]<<". CheckIsEdgePassThrowBlocked backward.\t after "<<TIMER_SEC(target)<<"s";
                    continue;
                }
                if(reached_target_phantomes.find(cur_node) != reached_target_phantomes.end())
                {
                    TIMER_STOP(target);
                    SimpleLogger().Write()<<"Bloking "<<i<<" "<<idxs[j]<<". reached_target_phantomes target.\t after "<<TIMER_SEC(target)<<"s";
                    continue;
                }
                nearest_forward_graph[i].push_back(idxs[j]);
                nearest_reverse_graph[idxs[j]].push_back(i);
                if(time_matrix(i, idxs[j]) > 0) reached_target_phantomes.insert(cur_node);
                TIMER_STOP(target);
                SimpleLogger().Write()<<"Nearest for "<<i<<" is "<<idxs[j]<<" found on "<<j<<" iteration.\t after "<<TIMER_SEC(target)<<"s";
                //if(nearest_graph[i].size() == NEAREST_RADIUS)
                //    threshold = 1 * time_matrix.at(i, idxs[j]);
            }
        }
    }
        
    /*void BuildNearestGraph(std::vector<std::set<PointID>> &nearest_graph, 
                           const WeightMatrix &time_matrix, 
                           const CrossNodesTable &cross_nodes_table, 
                           const HierarchyTree &forward_hierarchy_tree, 
                           const HierarchyTree &backward_hierarchy_tree,
                           const std::vector<std::set<NodeID>> &phantomes_for_point,
                           const std::set<NodeID> &all_phantomes,
                           const unsigned n) const
    {
        for(PointID i = 0; i < n; ++i)
        {
            TIMER_START(source);
            std::vector<PointID> idxs;
            for(PointID j = 0; j < n; ++j)
                if(i != j) idxs.push_back(j);
            std::sort(idxs.begin(), idxs.end(), [i, &time_matrix](PointID l, PointID r){
                return time_matrix.at(i, l)<time_matrix.at(i, r);
            });
            EdgeWeight threshold = std::numeric_limits<EdgeWeight>::max()-1;//idxs[n-2]/2;
            SimpleLogger().Write()<<"Initial threshold for "<<i<<"="<<threshold;
            std::set<NodeID> blocked_cross_nodes, 
                             all_blocked, 
                             reached_target_phantomes;
            for(unsigned j = 0; j < std::min(n - 1, 20u) && time_matrix.at(i, idxs[j]) <= threshold; ++j)
            {
                TIMER_START(target);
                NodeID cur_node = cross_nodes_table[i * n + idxs[j]];
                if(blocked_cross_nodes.find(cur_node) != blocked_cross_nodes.end())
                {
                    TIMER_STOP(target);
                    SimpleLogger().Write()<<"Bloking "<<i<<" "<<idxs[j]<<". blocked_cross_nodes.\t after "<<TIMER_SEC(target)<<"s"; 
                    continue;
                }
                                
                all_blocked.clear();
                std::set<NodeID> tempset;        
                std::set_union(phantomes_for_point[i].begin(), phantomes_for_point[i].end(),
                               phantomes_for_point[idxs[j]].begin(), phantomes_for_point[idxs[j]].end(),
                               std::inserter(tempset, tempset.begin()));
                std::set_difference(all_phantomes.begin(), all_phantomes.end(),
                                    tempset.begin(), tempset.end(),
                                    std::inserter(all_blocked, all_blocked.begin()));
                all_blocked.insert(blocked_cross_nodes.begin(), blocked_cross_nodes.end());
                bool is_first_node = true;
                
                SimpleLogger().Write(logDEBUG)<<"Test forward "<<i<<" "<<idxs[j];
                BOOST_ASSERT_MSG(forward_hierarchy_tree[i].find(cur_node) != forward_hierarchy_tree[i].end(), 
                                 ("forward_hierarchy_tree[" + boost::lexical_cast<std::string>(i) + "] has not key " + boost::lexical_cast<std::string>(cur_node)).c_str());
                while(cur_node != forward_hierarchy_tree[i].at(cur_node))
                {
                    SimpleLogger().Write(logDEBUG)<<cur_node;
                    if(CheckIsEdgePassThrowBlocked(forward_hierarchy_tree[i].at(cur_node), cur_node, all_blocked))
                    {
                        is_first_node=false;
                        break;
                    }
                    cur_node=forward_hierarchy_tree[i].at(cur_node);
                    BOOST_ASSERT_MSG(forward_hierarchy_tree[i].find(cur_node) != forward_hierarchy_tree[i].end(), 
                                 ("forward_hierarchy_tree[" + boost::lexical_cast<std::string>(i) + "] has not key " + boost::lexical_cast<std::string>(cur_node)).c_str());
                }
                if(!is_first_node)
                {
                    blocked_cross_nodes.insert(cross_nodes_table[i * n + idxs[j]]);
                    TIMER_STOP(target);
                    SimpleLogger().Write()<<"Bloking "<<i<<" "<<idxs[j]<<". CheckIsEdgePassThrowBlocked forward.\t after "<<TIMER_SEC(target)<<"s";
                    continue;
                }
                if(reached_target_phantomes.find(cur_node) != reached_target_phantomes.end())
                {
                    TIMER_STOP(target);
                    SimpleLogger().Write()<<"Bloking "<<i<<" "<<idxs[j]<<". reached_target_phantomes target.\t after "<<TIMER_SEC(target)<<"s";
                    continue;
                }
                
                cur_node = cross_nodes_table[i * n + idxs[j]];
                SimpleLogger().Write(logDEBUG)<<"Test backward "<<i<<" "<<idxs[j];
                BOOST_ASSERT_MSG(backward_hierarchy_tree[idxs[j]].find(cur_node) != backward_hierarchy_tree[idxs[j]].end(), 
                                 ("backward_hierarchy_tree[" + boost::lexical_cast<std::string>(i) + "] has not key " + boost::lexical_cast<std::string>(cur_node)).c_str());
                while(cur_node != backward_hierarchy_tree[idxs[j]].at(cur_node))
                {
                    SimpleLogger().Write(logDEBUG)<<cur_node;
                    if(CheckIsEdgePassThrowBlocked(cur_node, backward_hierarchy_tree[idxs[j]].at(cur_node), all_blocked))
                    {
                        is_first_node=false;
                        break;
                    }
                    cur_node=backward_hierarchy_tree[idxs[j]].at(cur_node);
                    BOOST_ASSERT_MSG(backward_hierarchy_tree[idxs[j]].find(cur_node) != backward_hierarchy_tree[idxs[j]].end(), 
                                     ("backward_hierarchy_tree[" + boost::lexical_cast<std::string>(i) + "] has not key " + boost::lexical_cast<std::string>(cur_node)).c_str());
                }
                if(!is_first_node) {
                    TIMER_STOP(target);
                    SimpleLogger().Write()<<"Bloking "<<i<<" "<<idxs[j]<<". CheckIsEdgePassThrowBlocked backward.\t after "<<TIMER_SEC(target)<<"s";
                    continue;
                }
                if(reached_target_phantomes.find(cur_node) != reached_target_phantomes.end())
                {
                    TIMER_STOP(target);
                    SimpleLogger().Write()<<"Bloking "<<i<<" "<<idxs[j]<<". reached_target_phantomes target.\t after "<<TIMER_SEC(target)<<"s";
                    continue;
                }
                nearest_graph[i].insert(idxs[j]);
                reached_target_phantomes.insert(cur_node);
                TIMER_STOP(target);
                SimpleLogger().Write()<<"Nearest for "<<i<<" is "<<idxs[j]<<" found on "<<j<<" iteration.\t after "<<TIMER_SEC(target)<<"s";
                if(nearest_graph[i].size() == NEAREST_RADIUS)
                    threshold = 1 * time_matrix.at(i, idxs[j]);
            }
            TIMER_STOP(source);
            SimpleLogger().Write()<<"Search nearest from "<<i<<" take "<<TIMER_SEC(source)<<"s";
        }
    }*/
    
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
    
    bool SetHasIntesection(const std::set<NodeID> &set1, const std::set<NodeID> &set2) const
    {
        auto first1=set1.begin(), 
             last1=set1.end(), 
             first2=set2.begin(), 
             last2=set2.end();
        while (first1!=last1 && first2!=last2)
        {
            if (*first1<*first2) ++first1;
            else if (*first2<*first1) ++first2;
            else return true;
        }
        return false;
    }
    void OutputChainsByStart (std::vector<char> &output,
                              const std::vector<std::list<PointID>> &nearest_graph) const
    {
        const unsigned n = nearest_graph.size();
        JSON::Object json_root;
        //JSON::Array chain_groups_array;
        JSON::Array chains_array;
        JSON::Array counts_array;
        JSON::Array nearest_array;
        for(PointID start=0; start<n; ++start)
        {
            /*for(const Chain &chain : chains_pull[start])
                if(chain.longest)
                {
                    JSON::Array chain_points_array;
                    chain_points_array.values.push_back(start);
                    for(const PointID point : chain.points_chain)
                        chain_points_array.values.push_back(point);
                    chains_array.values.push_back(chain_points_array);
                }
            JSON::Object point_counts_array;
            point_counts_array.values["all"]=points_use_count[0 * n + start];
            point_counts_array.values["short"]=points_use_count[1 * n + start];
            point_counts_array.values["longest"]=points_use_count[2 * n + start];
            counts_array.values.push_back(point_counts_array);*/
            JSON::Array nearest_array_row;
            for(const PointID point : nearest_graph[start])
                if(point>=0 && nearest_array_row.values.size()<2)
                    nearest_array_row.values.push_back(point);
            nearest_array.values.push_back(nearest_array_row);
        }
        json_root.values["chains"] = chains_array;
        json_root.values["counts"] = counts_array;
        json_root.values["nearest"] = nearest_array;
        json_root.values["n"] = n;
        JSON::render(output, json_root);
    }
    /*void RecursiveLookForChain(const PointID cur_point,
                               std::set<PointID> attendance_set, 
                               std::list<ChainRef> track_list,
                               std::vector<std::vector<Chain>> &chains_pull_ref,
                               std::vector<std::multimap<PointID, PointID>> &tails_list_ref,
                               std::vector<int> &points_use_count_ref,
                               const std::vector<std::set<PointID>> &nearest_graph,
                               const WeightMatrix &time_matrix,
                               const unsigned n) const
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
                chain_for_grow->Grow(cur_point, time_matrix, points_use_count_ref);
                for(PointID point : chain_from_cache.front()->points_chain)
                    chain_for_grow->Grow(point, time_matrix, points_use_count_ref);
            }
            return; //without recursion
        }
        else if(chain_from_cache.size() > 1)
        {
            for(ChainRef chain_for_grow : track_list)
                chain_for_grow->Grow(cur_point, time_matrix, points_use_count_ref);
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
            chain_for_grow->Grow(cur_point, time_matrix, points_use_count_ref);
        }
        
        attendance_set.insert(cur_point); // set current point as attended
        for(PointID point : attendance_set)
            ++points_use_count_ref[point];
        
        if(allowed_from_here.size() == 1)
        {
            // create new chain from current point
            chains_pull_ref[cur_point].emplace_back(cur_point, stopped_from_here, attendance_set.size()==1);
            ++points_use_count_ref[1 * n + cur_point];
            if(attendance_set.size()==1)
                ++points_use_count_ref[2 * n + cur_point];
            track_list.emplace_back(chains_pull_ref[cur_point].end() - 1); //track chain from current point
            RecursiveLookForChain(*allowed_from_here.begin(),
                                  attendance_set, 
                                  track_list,
                                  chains_pull_ref, 
                                  tails_list_ref, 
                                  points_use_count_ref, 
                                  nearest_graph, 
                                  time_matrix,
                                  n);
        }
        else if(allowed_from_here.size() > 1)
        {
            for(PointID next : allowed_from_here)
            {
                track_list.clear(); //stop traking
                // create new chain from current point
                chains_pull_ref[cur_point].emplace_back(cur_point, stopped_from_here, true);
                ++points_use_count_ref[1 * n + cur_point];
                track_list.emplace_back(chains_pull_ref[cur_point].end() - 1); //track chain from current point
                RecursiveLookForChain(next,
                                      attendance_set, 
                                      track_list,
                                      chains_pull_ref, 
                                      tails_list_ref, 
                                      points_use_count_ref, 
                                      nearest_graph, 
                                      time_matrix,
                                      n);
            }
        }
        // else size == 0 then return without recursion
    }
    
    void FilterChains(std::vector<std::vector<Chain>> &chains_pull)
    {
        std::set<Chain*, &Chain::CompareByEnd> ordered_by_end;
        for(PointID start=0; start<chains_pull.size(); ++start)
            for(Chain &chain : chains_pull[start])
                ordered_by_end.insert(&chain);
        auto prev_inter = std::partition_point(ordered_by_end.begin(), ordered_by_end.end(), [](Chain* c){ return !c->longest; });
        if(prev_inter != ordered_by_end.end()) {
            auto cur_inter = prev_inter + 1;
            for(; cur_inter!=ordered_by_end.end(); ++cur_inter)
                if(std::equal((*prev_inter)->points_chain.rbegin(), 
                              (*prev_inter)->points_chain.rend(),
                              (*cur_inter)->points_chain.rbegin()))
                    (*prev_inter)->longest = false;
        }
        for(PointID start=0; start<chains_pull.size(); ++start)
        {
            std::sort(chains_pull[start].begin(), chains_pull[start].end(), &Chain::CompareByStart);
            auto prev_inter = std::partition_point(chains_pull[start].begin(), chains_pull[start].end(), [](const Chain &c){ return !c.longest; });
            if(prev_inter != chains_pull[start].end()) {
                auto cur_inter = prev_inter + 1;
                for(; cur_inter!=chains_pull[start].end(); ++cur_inter)
                    if(std::equal(prev_inter->points_chain.begin(), 
                                  prev_inter->points_chain.end(),
                                  cur_inter->points_chain.begin()))
                        prev_inter->longest = false;
            }
        }
    }
    
    void OutputChainsByStart (std::vector<char> &output,
                              const std::vector<std::vector<Chain>> &chains_pull,
                              const std::vector<int> &points_use_count,
                              const std::vector<std::set<PointID>> &nearest_graph,
                              double time) const
    {
        const unsigned n = chains_pull.size();
        JSON::Object json_root;
        //JSON::Array chain_groups_array;
        JSON::Array chains_array;
        JSON::Array counts_array;
        JSON::Array nearest_array;
        for(PointID start=0; start<n; ++start)
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
            JSON::Object point_counts_array;
            point_counts_array.values["all"]=points_use_count[0 * n + start];
            point_counts_array.values["short"]=points_use_count[1 * n + start];
            point_counts_array.values["longest"]=points_use_count[2 * n + start];
            counts_array.values.push_back(point_counts_array);
            JSON::Array nearest_array_row;
            for(const PointID point : nearest_graph[start])
                nearest_array_row.values.push_back(point);
            nearest_array.values.push_back(nearest_array_row);
        }
        json_root.values["chains"] = chains_array;
        json_root.values["counts"] = counts_array;
        json_root.values["nearest"] = nearest_array;
        json_root.values["n"] = n;
        json_root.values["time"] = time;
        JSON::render(output, json_root);
    }
    
    bool CheckIsEdgePassThrowBlocked(NodeID from, NodeID to, std::set<NodeID> bloked) const
    {
        if(bloked.find(from) != bloked.end())
            return true;
        EdgeID edge = super::facade->FindEdgeInEitherDirection(from, to);
        BOOST_ASSERT(edge != SPECIAL_EDGEID);
        EdgeData ed = super::facade->GetEdgeData(edge);
        if(ed.shortcut)
            return CheckIsEdgePassThrowBlocked( from, ed.id, bloked)
                || CheckIsEdgePassThrowBlocked(ed.id,    to, bloked);
        else return false;
    }
    
    bool set_has_intesection(const std::set<NodeID> &set1, const std::set<NodeID> &set2) const
    {
        auto first1=set1.begin(), 
             last1=set1.end(), 
             first2=set2.begin(), 
             last2=set2.end();
        while (first1!=last1 && first2!=last2)
        {
            if (*first1<*first2) ++first1;
            else if (*first2<*first1) ++first2;
            else return true;
        }
        return false;
    }
    class Comp
    {
        const PointID from;
        const WeightMatrix &matrix; 
    public:
        Comp(const PointID from, const WeightMatrix &matrix)
            : from(from), matrix(matrix)
        {
        }
        bool operator()(PointID a, PointID b) const 
        {
            EdgeWeight a_dist, b_dist;
            if(a < 0) a_dist = matrix.at(1 - a, from);
            else a_dist = matrix.at(from, a);
            if(b < 0) b_dist = matrix.at(1 - b, from);
            else b_dist = matrix.at(from, b);
            return a_dist < b_dist;
        }
    };*/
    /*void BuildNearestGraph(//std::vector<std::set<PointID>> &nearest_graph, 
                           const ublas::matrix<EdgeWeight> &time_matrix, 
                           const ublas::matrix<EdgeWeight> &length_matrix, 
                           const ublas::matrix<NodeID> &cross_nodes_table, 
                           const HierarchyTree &forward_hierarchy_tree, 
                           const HierarchyTree &backward_hierarchy_tree,
                           const PhantomNodeArray &phantom_nodes_array) const
    {
        std::vector<std::set<PointID, Comp>> nearest_graph;
        const unsigned n = time_matrix.n;
        nearest_graph.reserve(n);
        for(PointID i = 0; i < n; ++i)
            nearest_graph.emplace_back(Comp(i, time_matrix));

        for(i = 0; i < n; ++i)
        {
            TIMER_START(source);
            std::vector<PointID> idxs;
            for(PointID j = 0; j < n; ++j)
                if(i != j) idxs.push_back(j);
            std::sort(idxs.begin(), idxs.end(), [i, &time_matrix](PointID l, PointID r){
                return time_matrix.at(i, l)<time_matrix.at(i, r);
            });
            EdgeWeight threshold = std::numeric_limits<EdgeWeight>::max()-1;//idxs[n-2]/2;
            SimpleLogger().Write()<<"Initial threshold for "<<i<<"="<<threshold;
            std::set<NodeID> blocked_cross_nodes, 
                             all_blocked, 
                             reached_target_phantomes;
            for(unsigned j = 0; j < std::min(n - 1, 20u) && time_matrix.at(i, j) <= threshold; ++j)
            {
                TIMER_START(target);
                NodeID cur_node = cross_nodes_table[i * n + j];
                //if(!cur_node)
                //{
                //    SimpleLogger().Write()<<"Bloking "<<i<<" "<<j<<". unreachable"; 
                //    continue;
                //}
                if(blocked_cross_nodes.find(cur_node) != blocked_cross_nodes.end())
                {
                    TIMER_STOP(target);
                    SimpleLogger().Write()<<"Bloking "<<i<<" "<<j<<". blocked_cross_nodes.\t after "<<TIMER_SEC(target)<<"s"; 
                    continue;
                }
                if(set_has_intesection(reached_target_phantomes, phantomes[j]))
                {
                    TIMER_STOP(target);
                    SimpleLogger().Write()<<"Bloking "<<i<<" "<<j<<". target phantom_node already reached.\t after "<<TIMER_SEC(target)<<"s";
                    continue;
                }
                                
                all_blocked.clear();
                std::set<NodeID> tempset;        
                std::set_union(phantomes[i].begin(), phantomes[i].end(),
                               phantomes[j].begin(), phantomes[j].end(),
                               std::inserter(tempset, tempset.begin()));
                std::set_difference(phantomes_all.begin(), phantomes_all.end(),
                                    tempset.begin(), tempset.end(),
                                    std::inserter(all_blocked, all_blocked.begin()));
                all_blocked.insert(blocked_cross_nodes.begin(), blocked_cross_nodes.end());
                bool is_first_node = true;
                
                SimpleLogger().Write(logDEBUG)<<"Test forward "<<i<<" "<<j;
                BOOST_ASSERT_MSG(forward_hierarchy_tree[i].find(cur_node) != forward_hierarchy_tree[i].end(), 
                                 ("forward_hierarchy_tree[" + boost::lexical_cast<std::string>(i) + "] has not key " + boost::lexical_cast<std::string>(cur_node)).c_str());
                while(cur_node != forward_hierarchy_tree[i].at(cur_node))
                {
                    SimpleLogger().Write(logDEBUG)<<cur_node;
                    if(CheckIsEdgePassThrowBlocked(forward_hierarchy_tree[i].at(cur_node), cur_node, all_blocked))
                    {
                        is_first_node=false;
                        break;
                    }
                    cur_node=forward_hierarchy_tree[i].at(cur_node);
                    BOOST_ASSERT_MSG(forward_hierarchy_tree[i].find(cur_node) != forward_hierarchy_tree[i].end(), 
                                 ("forward_hierarchy_tree[" + boost::lexical_cast<std::string>(i) + "] has not key " + boost::lexical_cast<std::string>(cur_node)).c_str());
                }
                if(!is_first_node)
                {
                    blocked_cross_nodes.insert(cross_nodes_table[i * n + j]);
                    TIMER_STOP(target);
                    SimpleLogger().Write()<<"Bloking "<<i<<" "<<j<<". CheckIsEdgePassThrowBlocked forward.\t after "<<TIMER_SEC(target)<<"s";
                    continue;
                }
                if(reached_target_phantomes.find(cur_node) != reached_target_phantomes.end())
                {
                    TIMER_STOP(target);
                    SimpleLogger().Write()<<"Bloking "<<i<<" "<<j<<". reached_target_phantomes target.\t after "<<TIMER_SEC(target)<<"s";
                    continue;
                }
                
                cur_node = cross_nodes_table[i * n + j];
                SimpleLogger().Write(logDEBUG)<<"Test backward "<<i<<" "<<j;
                BOOST_ASSERT_MSG(backward_hierarchy_tree[j].find(cur_node) != backward_hierarchy_tree[j].end(), 
                                 ("backward_hierarchy_tree[" + boost::lexical_cast<std::string>(i) + "] has not key " + boost::lexical_cast<std::string>(cur_node)).c_str());
                while(cur_node != backward_hierarchy_tree[j].at(cur_node))
                {
                    SimpleLogger().Write(logDEBUG)<<cur_node;
                    if(CheckIsEdgePassThrowBlocked(cur_node, backward_hierarchy_tree[j].at(cur_node), all_blocked))
                    {
                        is_first_node=false;
                        break;
                    }
                    cur_node=backward_hierarchy_tree[j].at(cur_node);
                    BOOST_ASSERT_MSG(backward_hierarchy_tree[j].find(cur_node) != backward_hierarchy_tree[j].end(), 
                                     ("backward_hierarchy_tree[" + boost::lexical_cast<std::string>(i) + "] has not key " + boost::lexical_cast<std::string>(cur_node)).c_str());
                }
                if(!is_first_node) {
                    TIMER_STOP(target);
                    SimpleLogger().Write()<<"Bloking "<<i<<" "<<j<<". CheckIsEdgePassThrowBlocked backward.\t after "<<TIMER_SEC(target)<<"s";
                    continue;
                }
                if(reached_target_phantomes.find(cur_node) != reached_target_phantomes.end())
                {
                    TIMER_STOP(target);
                    SimpleLogger().Write()<<"Bloking "<<i<<" "<<j<<". reached_target_phantomes target.\t after "<<TIMER_SEC(target)<<"s";
                    continue;
                }
                nearest_graph[i].insert(j);
                reached_target_phantomes.insert(cur_node);
                TIMER_STOP(target);
                SimpleLogger().Write()<<"Nearest for "<<i<<" is "<<j<<" found on "<<j<<" iteration.\t after "<<TIMER_SEC(target)<<"s";
                if(nearest_graph[i].size() == NEAREST_RADIUS)
                    threshold = 1 * time_matrix.at(i, j);
            }
            TIMER_STOP(source);
            SimpleLogger().Write()<<"Search nearest from "<<i<<" take "<<TIMER_SEC(source)<<"s";
        }
    }*/

    void ForwardRoutingStep(const unsigned source_id,
                            const unsigned n,
                            QueryHeap &query_heap,
                            LengthMap &length_map,
                            const SearchSpaceWithBuckets &search_space_with_buckets,
                            HierarchyTree &hierarchy_tree,
                            const HierarchyTree &backward_hierarchy_tree,
                            ublas::matrix<EdgeWeight> &time_matrix, 
                            ublas::matrix<EdgeWeight> &length_matrix, 
                            ublas::matrix<NodeID> &cross_nodes_table, 
                            const TransportRestriction &tr) const
    {
        const NodeID node = query_heap.DeleteMin();
        const int source_distance = query_heap.GetKey(node);

        hierarchy_tree[source_id][node] = query_heap.GetData(node).parent;

        
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
                const EdgeWeight new_distance = source_distance + target_distance;
                const EdgeWeight current_distance = time_matrix(source_id, target_id);
                if(new_distance >= 0 && new_distance < current_distance)
                {
                    cross_nodes_table(source_id, target_id) = node;
                    time_matrix(source_id, target_id) = new_distance;
                    length_matrix(source_id, target_id) = abs(length_map[node] + current_bucket.length);
                }
            }
        }
        if (StallAtNode<true>(node, source_distance, query_heap, tr))
        {
            return;
        }
        RelaxOutgoingEdges<true>(node, source_distance, query_heap, length_map, tr);
    }

    void BackwardRoutingStep(const unsigned target_id,
                             QueryHeap &query_heap,
                             LengthMap &length_map,
                             SearchSpaceWithBuckets &search_space_with_buckets, 
                             HierarchyTree &hierarchy_tree,
                             const TransportRestriction &tr) const
    {
        const NodeID node = query_heap.DeleteMin();
        const int target_distance = query_heap.GetKey(node);

        // store settled nodes in search space bucket
        search_space_with_buckets[node].emplace_back(target_id, target_distance, length_map[node], query_heap.GetData(node).parent);
        hierarchy_tree[target_id][node] = query_heap.GetData(node).parent;

        if (StallAtNode<false>(node, target_distance, query_heap, tr))
        {
            return;
        }

        RelaxOutgoingEdges<false>(node, target_distance, query_heap, length_map, tr);
    }

    template <bool forward_direction>
    inline void
    RelaxOutgoingEdges(const NodeID node, const EdgeWeight distance, QueryHeap &query_heap, LengthMap &length_map, const TransportRestriction &tr) const
    {
        for (auto edge : super::facade->GetAdjacentEdgeRange(node))
        {
            const auto &data = super::facade->GetEdgeData(edge);
            if(tr.IsEdgeRestricted(data)) continue;
            const bool direction_flag = (forward_direction ? data.forward : data.backward);
            if(direction_flag)
            {
                const NodeID to = super::facade->GetTarget(edge);
                const int edge_weight = data.distance;

                BOOST_ASSERT_MSG(edge_weight > 0, "edge_weight invalid");
                const int to_distance = distance + edge_weight;

                // New Node discovered -> Add to Heap + Node Info Storage
                if (!query_heap.WasInserted(to))
                {
                    query_heap.Insert(to, to_distance, node);
                    length_map[to] = length_map[node] + data.length;
                }
                // Found a shorter Path -> Update distance
                else if (to_distance < query_heap.GetKey(to))
                {
                    // new parent
                    query_heap.GetData(to).parent = node;
                    query_heap.DecreaseKey(to, to_distance);
                    length_map[to] = length_map[node] + data.length;
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
