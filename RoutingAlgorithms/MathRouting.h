#ifndef MATH_ROUTING_H
#define MATH_ROUTING_H

#include "BasicRoutingInterface.h"
#include "../DataStructures/JSONContainer.h"
#include "../DataStructures/SearchEngineData.h"
#include "../typedefs.h"
#include "../Util/TimingUtil.h"

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
        PointID target_id; // essentially a row in the distance matrix
        EdgeWeight distance;
        NodeID parent;
        NodeBucket(const unsigned target_id, const EdgeWeight distance, const NodeID parent)
            : target_id(target_id), distance(distance), parent(parent)
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
        void Grow(PointID point, const WeightMatrix &weight_matrix, std::vector<int> &points_use_count)
        {
            forward_distance += weight_matrix.at(GetLast(), point);
            reverse_distance += weight_matrix.at(point, GetLast());
            points_chain.push_back(point);
            points_set.insert(point);
            ++points_use_count[1 * weight_matrix.number_of_locations + point];
            if(longest)
                ++points_use_count[2 * weight_matrix.number_of_locations + point];
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
        
        static bool CompareByStart(const Chain &left, const Chain &right)
        {
            return !left.longest || 
                std::lexicographical_compare( left.points_chain.begin(),  left.points_chain.end(), 
                                             right.points_chain.begin(), right.points_chain.end());
        }
        
        static bool CompareByEnd(const Chain *left, const Chain *right)
        {
            if(!left->longest) return true;
            bool res = std::lexicographical_compare( left->points_chain.rbegin(),  left->points_chain.rend(), 
                                                    right->points_chain.rbegin(), right->points_chain.rend());
            
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
    typedef std::vector<std::unordered_map<NodeID, NodeID>> HierarchyTree;
    typedef std::vector<NodeID>  CrossNodesTable;
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
        TIMER_START(process);
        const unsigned number_of_locations=phantom_nodes_array.size();
        std::shared_ptr<std::vector<unsigned>> result_table =
            std::make_shared<std::vector<unsigned>>(number_of_locations,
                    std::numeric_limits<unsigned>::max());
        WeightMatrix time_matrix(number_of_locations);

        engine_working_data.InitializeOrClearFirstThreadLocalStorage(
            super::facade->GetNumberOfNodes());

        QueryHeap &query_heap = *(engine_working_data.forwardHeap);

        SearchSpaceWithBuckets backward_search_space_with_buckets;
        HierarchyTree backward_hierarchy_tree(number_of_locations), forward_hierarchy_tree(number_of_locations);
        CrossNodesTable cross_nodes_table(number_of_locations*number_of_locations);

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
                BackwardRoutingStep(target_id, query_heap, backward_search_space_with_buckets, backward_hierarchy_tree, tr);
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
                                   forward_hierarchy_tree,
                                   backward_hierarchy_tree,
                                   cross_nodes_table,
                                   time_matrix,
                                   tr);
            }

            ++source_id;
        }
        BOOST_ASSERT(source_id == target_id);
        
        // 0 - all
        // 1 - short
        // 3 - longest
        std::vector<int> points_use_count(3 * number_of_locations);
        
        std::vector<std::set<PointID>> nearest_graph(number_of_locations);
        
        std::vector<std::vector<Chain>> chains_pull(number_of_locations);
        
        //i=end>>>|       one-before-end|  |start|
        std::vector<std::multimap<PointID, PointID>> tails_list(number_of_locations);
        
        // Building Nearest Graph
        BuildNearestGraph(nearest_graph, 
                          time_matrix, 
                          cross_nodes_table, 
                          forward_hierarchy_tree, 
                          backward_hierarchy_tree,
                          phantom_nodes_array);
        
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
        //FilterChains(chains_pull);
        TIMER_STOP(process);
        OutputChainsByStart(output, chains_pull, points_use_count, nearest_graph, TIMER_SEC(process));
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
            ++points_use_count_ref[1 * number_of_locations + cur_point];
            if(attendance_set.size()==1)
                ++points_use_count_ref[2 * number_of_locations + cur_point];
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
                ++points_use_count_ref[1 * number_of_locations + cur_point];
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
    
    /*void FilterChains(std::vector<std::vector<Chain>> &chains_pull)
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
    }*/
    
    void OutputChainsByStart (std::vector<char> &output,
                              const std::vector<std::vector<Chain>> &chains_pull,
                              const std::vector<int> &points_use_count,
                              const std::vector<std::set<PointID>> &nearest_graph,
                              double time) const
    {
        const unsigned number_of_locations = chains_pull.size();
        JSON::Object json_root;
        //JSON::Array chain_groups_array;
        JSON::Array chains_array;
        JSON::Array counts_array;
        JSON::Array nearest_array;
        for(PointID start=0; start<number_of_locations; ++start)
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
            point_counts_array.values["all"]=points_use_count[0 * number_of_locations + start];
            point_counts_array.values["short"]=points_use_count[1 * number_of_locations + start];
            point_counts_array.values["longest"]=points_use_count[2 * number_of_locations + start];
            counts_array.values.push_back(point_counts_array);
            JSON::Array nearest_array_row;
            for(const PointID point : nearest_graph[start])
                nearest_array_row.values.push_back(point);
            nearest_array.values.push_back(nearest_array_row);
        }
        json_root.values["chains"] = chains_array;
        json_root.values["counts"] = counts_array;
        json_root.values["nearest"] = nearest_array;
        json_root.values["n"] = number_of_locations;
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
    
    void BuildNearestGraph(std::vector<std::set<PointID>> &nearest_graph, 
                           const WeightMatrix &time_matrix, 
                           const CrossNodesTable &cross_nodes_table, 
                           const HierarchyTree &forward_hierarchy_tree, 
                           const HierarchyTree &backward_hierarchy_tree,
                           const PhantomNodeArray &phantom_nodes_array) const
    {
        const unsigned n = time_matrix.number_of_locations;
        std::vector<std::set<NodeID>> phantomes(n);
        std::set<NodeID> phantomes_all;
        PointID i = 0;
        for (const std::vector<PhantomNode> &phantom_node_vector : phantom_nodes_array)
        {
            for (const PhantomNode &phantom_node : phantom_node_vector)
            {
                if (SPECIAL_NODEID != phantom_node.forward_node_id)
                {
                    phantomes[i].insert(phantom_node.forward_node_id);
                    phantomes_all.insert(phantom_node.forward_node_id);
                    SimpleLogger().Write(logDEBUG) << "Phantome for "<<i<<" is "<<phantom_node.forward_node_id;
                }
                if (SPECIAL_NODEID != phantom_node.reverse_node_id)
                {
                    phantomes[i].insert(phantom_node.reverse_node_id);
                    phantomes_all.insert(phantom_node.reverse_node_id);
                    SimpleLogger().Write(logDEBUG) << "Phantome for "<<i<<" is "<<phantom_node.reverse_node_id;
                }
            }
            ++i;
        }
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
                NodeID cur_node = cross_nodes_table[i * n + idxs[j]];
                /*if(!cur_node)
                {
                    SimpleLogger().Write()<<"Bloking "<<i<<" "<<idxs[j]<<". unreachable"; 
                    continue;
                }*/
                if(blocked_cross_nodes.find(cur_node) != blocked_cross_nodes.end())
                {
                    TIMER_STOP(target);
                    SimpleLogger().Write()<<"Bloking "<<i<<" "<<idxs[j]<<". blocked_cross_nodes.\t after "<<TIMER_SEC(target)<<"s"; 
                    continue;
                }
                if(set_has_intesection(reached_target_phantomes, phantomes[idxs[j]]))
                {
                    TIMER_STOP(target);
                    SimpleLogger().Write()<<"Bloking "<<i<<" "<<idxs[j]<<". target phantom_node already reached.\t after "<<TIMER_SEC(target)<<"s";
                    continue;
                }
                                
                all_blocked.clear();
                std::set<NodeID> tempset;        
                std::set_union(phantomes[i].begin(), phantomes[i].end(),
                               phantomes[idxs[j]].begin(), phantomes[idxs[j]].end(),
                               std::inserter(tempset, tempset.begin()));
                std::set_difference(phantomes_all.begin(), phantomes_all.end(),
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
                    threshold = 1.1 * time_matrix.at(i, idxs[j]);
            }
            TIMER_STOP(source);
            SimpleLogger().Write()<<"Search nearest from "<<i<<" take "<<TIMER_SEC(source)<<"s";
                
        }
    }

    void ForwardRoutingStep(const unsigned source_id,
                            const unsigned number_of_locations,
                            QueryHeap &query_heap,
                            const SearchSpaceWithBuckets &backward_search_space_with_buckets,
                            HierarchyTree &hierarchy_tree,
                            const HierarchyTree &backward_hierarchy_tree,
                            CrossNodesTable &cross_nodes_table,
                            WeightMatrix &time_matrix,
                            const TransportRestriction &tr) const
    {
        const NodeID node = query_heap.DeleteMin();
        const int source_distance = query_heap.GetKey(node);
        //const int source_length = length_map[node];

        hierarchy_tree[source_id][node] = query_heap.GetData(node).parent;
        BOOST_ASSERT_MSG(hierarchy_tree[source_id].find(query_heap.GetData(node).parent) != hierarchy_tree[source_id].end(), 
                         ("hierarchy_tree[" + boost::lexical_cast<std::string>(source_id) + "] "
                          "has not parent " + boost::lexical_cast<std::string>(query_heap.GetData(node).parent) + 
                          "for node " + boost::lexical_cast<std::string>(node)).c_str());

        
        // check if each encountered node has an entry
        const auto bucket_iterator = backward_search_space_with_buckets.find(node);
        // iterate bucket if there exists one
        if (bucket_iterator != backward_search_space_with_buckets.end())
        {
            const std::vector<NodeBucket> &bucket_list = bucket_iterator->second;
            for (const NodeBucket &current_bucket : bucket_list)
            {
                // get target id from bucket entry
                const unsigned target_id = current_bucket.target_id;
                const int target_distance = current_bucket.distance;
                const EdgeWeight new_distance = source_distance + target_distance;
                const EdgeWeight current_distance = time_matrix.at(source_id, target_id);
                if(new_distance >= 0 && new_distance < current_distance)
                {
                    cross_nodes_table[source_id * number_of_locations + target_id] = node;
                    BOOST_ASSERT_MSG(backward_hierarchy_tree[target_id].find(node) != backward_hierarchy_tree[source_id].end(), 
                                    ("backward_hierarchy_tree[" + boost::lexical_cast<std::string>(target_id) + "] "
                                    "has not node " + boost::lexical_cast<std::string>(node)).c_str());
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
                             HierarchyTree &hierarchy_tree,
                             const TransportRestriction &tr) const
    {
        const NodeID node = query_heap.DeleteMin();
        const int target_distance = query_heap.GetKey(node);

        // store settled nodes in search space bucket
        backward_search_space_with_buckets[node].emplace_back(target_id, target_distance, query_heap.GetData(node).parent);
        hierarchy_tree[target_id][node] = query_heap.GetData(node).parent;
        BOOST_ASSERT_MSG(hierarchy_tree[target_id].find(query_heap.GetData(node).parent) != hierarchy_tree[target_id].end(), 
                         ("hierarchy_tree[" + boost::lexical_cast<std::string>(target_id) + "] "
                          "has not parent " + boost::lexical_cast<std::string>(query_heap.GetData(node).parent) + 
                          "for node " + boost::lexical_cast<std::string>(node)).c_str());

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
