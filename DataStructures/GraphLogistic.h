#ifndef GRAPG_LOGISTIC_H
#define GRAPG_LOGISTIC_H

#include <boost/assert.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
//#include <boost/iterator/iterator_concepts.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include "JSONContainer.h"

#include <vector>
#include <list>
#include <algorithm>
#include <set>
#include <map>
#include <utility>
#include <ostream>
#include <memory>

namespace ublas = boost::numeric::ublas;

template<typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
    std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, " "));
    return out;
}

class GraphLogistic
{
    typedef unsigned PointID;
    typedef unsigned CoreID;
    typedef unsigned ChainID;
    typedef std::pair<double, double> Coordinate;
    typedef std::vector<PointID>::iterator TourBound;
    typedef std::shared_ptr<std::pair<TourBound, TourBound>> TourBoundsPair;
    typedef enum ChainType { POINT, MDIR, BDIR } ChainType;
    typedef std::pair<ChainType, std::list<PointID>> Chain;
    
    const int UNREACHED_WEIGHT = std::numeric_limits<int>::max();
    const unsigned NEAREST_RADIUS = 2;
    const unsigned MIN_CORE_RADIUS = 3;
    const unsigned MOVING_MEAN_RADIUS = 0;
    const float MOVING_MEAN_FACTOR = 1.9;
    const unsigned MAX_TOUR_LENGTH = 20;
    
    unsigned n_points, n_cores, n_chains;
    ublas::matrix<int> time_matrix, length_matrix;
    std::vector<Coordinate> coordinates;
    std::vector<std::set<PointID>> full_forward_graph, full_reverse_graph; 
    std::vector<std::set<PointID>> nearest_forward_graph, nearest_reverse_graph;
    std::vector<std::set<ChainID>> chain_forward_graph, chain_reverse_graph;
    std::vector<bool> tails_dfs_visited;
    std::vector<bool> tarjan_dfs_visited;
    std::list<PointID> tarjan_stack;
    std::vector<int> tarjan_lowlink;
    int tarjan_time;
    std::vector<std::set<PointID>> core_points, core_tails;
    std::vector<int> core_distances;
    std::vector<PointID> core_start_points;
    std::map<PointID, CoreID> point_cores;
    std::vector<std::set<CoreID>> tail_cores;
    std::vector<std::set<PointID>> tail_gates;
    std::map<ChainID, Chain> chains_index;
    std::list<std::vector<PointID>> tours;
    std::set<PointID> attended_points;
    std::set<CoreID> attended_cores;
    unsigned begin_offset, end_offset;
    
    void findCorePoints() 
    {
        tarjan_time = 0;
        tarjan_stack.clear();
        tarjan_dfs_visited.assign(n_points, false);
        tarjan_lowlink.assign(n_points, 0);
        tail_cores.assign(n_points, std::set<CoreID>());
        tail_gates.assign(n_points, std::set<PointID>());
        core_points.clear();
        core_points.reserve(n_points * 0.6);

        for (PointID u = 0; u < n_points; u++)
            if (!tarjan_dfs_visited[u])
                tarjanDFS(u);
        
        n_cores = core_points.size();
    }
    
    void tarjanDFS(PointID u) 
    {
        tarjan_lowlink[u] = tarjan_time++;
        tarjan_dfs_visited[u] = true;
        tarjan_stack.push_back(u);
        bool isComponentRoot = true;

        for (PointID v : nearest_forward_graph[u]) {
            if (!tarjan_dfs_visited[v])
                tarjanDFS(v);
            if (tarjan_lowlink[u] > tarjan_lowlink[v]) {
                tarjan_lowlink[u] = tarjan_lowlink[v];
                isComponentRoot = false;
            }
        }

        if (isComponentRoot) {
            std::set<PointID> core;
            while (true) {
                PointID x = tarjan_stack.back();
                tarjan_stack.pop_back();
                if(x != 0) core.insert(x);
                tarjan_lowlink[x] = std::numeric_limits<int>::max();
                if (x == u)
                    break;
            }
            if(core.size() >= MIN_CORE_RADIUS)
            {
                core_points.push_back(core);
                const CoreID cur_core = core_points.size() - 1;
                int min_dist = UNREACHED_WEIGHT;
                PointID min;
                for(PointID p : core)
                {
                    point_cores[p] = cur_core;
                    tail_cores[p].insert(cur_core);
                    tail_gates[p].insert(p);
                    //tails_dfs_visited[p] = true;
                    if(min_dist > time_matrix(0, p))
                    {
                        min_dist = time_matrix(0, p);
                        min = p;
                    }
                }
                core_distances.push_back(min_dist);
                core_start_points.push_back(min);
                //cores_order.insert(cur_core);
            }
        }
    }
    
    void findCoreForwardTails()
    {
        tails_dfs_visited.assign(n_points, false);
        core_tails.assign(n_cores, std::set<PointID>());

        for (PointID u = 0; u < n_points; u++)
            if (!tails_dfs_visited[u])
                forwardTailsDFS(u);
        
        rebindTail(boost::make_counting_iterator(0u), 
                   boost::make_counting_iterator(n_points),
                   [&](PointID p) { return tail_cores[p].empty(); });
            
    }

    void forwardTailsDFS(PointID u) 
    {
        tails_dfs_visited[u] = true;
        for (PointID v : nearest_forward_graph[u]) {
            if (!tails_dfs_visited[v] && point_cores.find(v) == point_cores.end())
                forwardTailsDFS(v);
            for(const CoreID &c : tail_cores[v])
            {
                core_tails[c].insert(u);
                tail_cores[u].insert(c);
            }
            if(point_cores.find(v) == point_cores.end() || point_cores.find(u) == point_cores.end())
                tail_gates[u].insert(tail_gates[v].begin(), tail_gates[v].end());
        }
    }
    
    template<class Iter, class Filter>
    void rebindTail(Iter begin, Iter end, Filter filter)
    {
        for(Iter p = begin; p != end; p++)
            if(filter(*p))
            {
                //for(CoreID c : tail_cores[*p])
                //    core_tails[c].erase(*p);
                //tail_cores[*p].clear();
                int min_w = UNREACHED_WEIGHT;
                PointID min;
                for(const PointID i : full_forward_graph[*p])
                    if(!filter(i) && min_w > time_matrix(*p, i))
                    {
                        min_w = time_matrix(*p, i);
                        min = i;
                    }
                if(min_w != UNREACHED_WEIGHT)
                {
                    for(const CoreID &c : tail_cores[min])
                    {
                        core_tails[c].insert(*p);
                        tail_cores[*p].insert(c);
                    }
                    tail_gates[*p].insert(tail_gates[min].begin(), tail_gates[min].end());
                }
            }
    }
    
    void findCoreBackwardTails()
    {
        tails_dfs_visited.assign(n_points, false);
        for (CoreID c = 0; c < n_cores; c++)
            for(const PointID g : core_points[c])
                for(const PointID v : nearest_forward_graph[g])
                    if(!tails_dfs_visited[v] && core_points[c].find(v) == core_points[c].end())
                        backwardTailsDFS(c, g, v);
    }

    void backwardTailsDFS(const CoreID c, const PointID g, const PointID u) 
    {
        tails_dfs_visited[u] = true;
        core_tails[c].insert(u);
        tail_cores[u].insert(c);
        tail_gates[u].insert(g);
        
        if (point_cores.find(u) == point_cores.end())
            for (const PointID v : nearest_forward_graph[u])
                if (!tails_dfs_visited[v] && point_cores.find(v) == point_cores.end())
                    backwardTailsDFS(c, g, v);
    }

    void linkUnattachedCores()
    {
        if(n_cores > 1)
            for(CoreID c = 0; c < n_cores; ++c)
            {
                bool unattached = true;
                for(PointID t : core_tails[c])
                    if(tail_cores[t].size() > 1) {
                        unattached = false;
                        break;
                    }
                if(unattached)
                {
                    PointID min_from, min_to;
                    int min_weight = UNREACHED_WEIGHT;
                    for(PointID from : core_tails[c])
                        for(PointID to : full_forward_graph[from])
                            if(tail_cores[to].find(c) == tail_cores[to].end() && min_weight > time_matrix(from, to))
                            {
                                min_weight = time_matrix(from, to);
                                min_from = from;
                                min_to = to;
                            }
                    tail_cores[min_to].insert(c);
                    tail_gates[min_to].insert(tail_gates[min_from].begin(), tail_gates[min_from].end());
                    tail_gates[min_from].insert(tail_gates[min_to].begin(), tail_gates[min_to].end());
                }
            }
    }
    
    float vectprod(std::pair<double,double> &A, std::pair<double,double> &B, std::pair<double,double> &C)
    { 
        return (B.first-A.first) * (C.second-B.second) - (B.second-A.second) * (C.first-B.first); 
    }
    
    void coreCWVisitation(std::vector<PointID> &tour, CoreID cur_core, PointID from, PointID first)
    {
        PointID second;
        PointID cur_node = first;
        PointID prev_node = from;
        attended_points.insert(cur_node);
        std::set<std::pair<PointID, PointID>> attended_edges;
        tour.push_back(cur_node);
        while(true)
        {
            auto i = full_forward_graph[cur_node].begin();
            while(core_points[cur_core].find(*i) == core_points[cur_core].end()) ++i;
            PointID next=*i;
            bool leftmost=vectprod(coordinates[prev_node], coordinates[cur_node], coordinates[next])<0;
            for(++i; i!=full_forward_graph[cur_node].end(); ++i)
            {
                if(core_points[cur_core].find(*i) == core_points[cur_core].end()) continue;
                if(!leftmost && vectprod(coordinates[prev_node], coordinates[cur_node], coordinates[*i])<0){
                    next=*i;
                    leftmost=true;
                }
                else if(*i!=prev_node && vectprod(coordinates[cur_node], coordinates[next], coordinates[*i])<0
                    && (!leftmost || vectprod(coordinates[prev_node], coordinates[cur_node], coordinates[*i])<0))
                        next=*i;
            }
            prev_node=cur_node;
            cur_node=next;
            auto cur_edge = std::make_pair(prev_node, cur_node);
            if(attended_edges.find(cur_edge) == attended_edges.end())
                attended_edges.insert(cur_edge);
            else break;
            if(attended_points.find(next) == attended_points.end())
            {
                attended_points.insert(cur_node);
                tour.push_back(cur_node);
            }
            //if(tours.back().size() >= MAX_TOUR_LENGTH)
            //    tours.emplace_back();
        }
    }
    
    int cost(const std::vector<PointID> &tour)
    {
        int sum=time_matrix(0, tour.front()) + time_matrix(tour.back(), 0);
        for(unsigned i=1; i<tour.size(); ++i)
        {
            sum+=time_matrix(tour[i - 1], tour[i]);
            for(unsigned j=1; j<=MOVING_MEAN_RADIUS; ++j)
            {
                if(int(i - j) > 0)
                    sum += time_matrix(tour[i - j - 1], tour[i]) / (j * MOVING_MEAN_FACTOR);
                if(i + j < tour.size())
                    sum += time_matrix(tour[i - 1], tour[i + j]) / (j * MOVING_MEAN_FACTOR);
            }
        }
        return sum;
    }
    
    bool checkOpt(const std::vector<PointID> &tempTour, std::vector<PointID> &bestTour, int &bestCost)
    {
        int tempCost = cost(tempTour);
        if(bestCost > tempCost)
        {
            std::cout << "route " << tempTour << " cost " << tempCost << " is best!" << std::endl;
            bestTour = std::move(tempTour);
            bestCost = tempCost;
            return true;
        }
        return false;
    }
    
    template<class Action> 
    void do1opt(std::vector<PointID> &tour, Action action, TourBound lbound, TourBound rbound, std::vector<PointID> &bestTour, int &bestCost, bool &forward)
    {
        std::cout << "opt1" << std::endl; 
        for(auto iter=lbound; iter <= rbound; ++iter)
        {
            
            std::vector<PointID> tempTour;
            tempTour.reserve(tour.size());
            
            std::copy(tour.begin(), iter, std::back_inserter(tempTour));
            action(tempTour);
            std::copy(iter, tour.end(), std::back_inserter(tempTour));
            
            
            if(checkOpt(tempTour, bestTour, bestCost))
                forward = true;
        }
    }
    
    template<class Action> 
    void do2opt(std::vector<PointID> &tour, Action action, TourBound lbound, TourBound rbound, std::vector<PointID> &bestTour, int &bestCost, bool &forward)
    {
        std::cout << "opt2" << std::endl;
        for(auto iter1=lbound; iter1 <= rbound; ++iter1)
            for(auto iter2=lbound; iter2 <= rbound; ++iter2)
            {
                std::vector<PointID> tempTour;
                tempTour.reserve(tour.size());
                
                if(iter2 > iter1 + 1)
                {
                    std::copy(tour.begin(), iter1, std::back_inserter(tempTour));
                    action(tempTour);
                    std::reverse_copy(iter1, iter2, std::back_inserter(tempTour));
                    std::copy(iter2, tour.end(), std::back_inserter(tempTour));
                }
                else if(iter2 < iter1)
                {
                    std::reverse_copy(iter1, tour.end(), std::back_inserter(tempTour));
                    std::copy(iter2, iter1, std::back_inserter(tempTour));
                    action(tempTour);
                    std::reverse_copy(tour.begin(), iter2, std::back_inserter(tempTour));
                }
                else continue;
                
                if(checkOpt(tempTour, bestTour, bestCost))
                    forward = iter2 > iter1 + 1;
            }
    }
    
    template<class Action> 
    void do3opt(std::vector<PointID> &tour, Action action, TourBound lbound, TourBound rbound, std::vector<PointID> &bestTour, int &bestCost, bool &forward)
    {
        std::cout << "opt3" << std::endl;
        for(auto iter1=lbound; iter1 <= rbound; ++iter1)
            for(auto iter2=iter1 + 2; iter2 <= rbound; ++iter2)
                for(auto iter3=lbound; iter3 <= rbound; ++iter3)
                {
                    std::vector<PointID> tempTour;
                    tempTour.reserve(tour.size());
                    
                    if(iter3 > iter2 + 1)
                    {
                        std::copy(tour.begin(), iter1, std::back_inserter(tempTour));
                        action(tempTour);
                        std::reverse_copy(iter1, iter2, std::back_inserter(tempTour));
                        std::reverse_copy(iter2, iter3, std::back_inserter(tempTour));
                        std::copy(iter3, tour.end(), std::back_inserter(tempTour));
                        
                    }
                    else if(iter3 < iter1)
                    {
                        std::reverse_copy(iter2, tour.end(), std::back_inserter(tempTour));
                        std::copy(iter3, iter1, std::back_inserter(tempTour));
                        action(tempTour);
                        std::reverse_copy(iter1, iter2, std::back_inserter(tempTour));
                        std::reverse_copy(tour.begin(), iter3, std::back_inserter(tempTour));
                    }
                    else continue;
                    
                    if(checkOpt(tempTour, bestTour, bestCost))
                        forward = iter3 > iter2 + 1;
                }
    }
    
    template<class Action> 
    bool doAllOpts(std::vector<PointID> &tour, Action action, TourBound lbound, TourBound rbound, std::vector<PointID> &bestTour, int &bestCost, bool &forward)
    {
        const int prevCost = bestCost;
        do1opt(tour, action, lbound, rbound, bestTour, bestCost, forward);
        do2opt(tour, action, lbound, rbound, bestTour, bestCost, forward);
        do3opt(tour, action, lbound, rbound, bestTour, bestCost, forward);
        return bestCost < prevCost;
    }

    bool tryInsert(PointID point, TourBound lbound, TourBound rbound)
    {
        std::cout << "try insert " << point << " somewhere between " << *lbound << " and " << *(rbound - 1) << std::endl;
        std::vector<PointID> bestTour;
        int bestCost = UNREACHED_WEIGHT;
        bool forward;
        
        auto insertPoint = [point](std::vector<PointID> &tour){ tour.push_back(point); };
        auto doNothing = [](std::vector<PointID> &tour){ };
        bool success = doAllOpts(tours.back(), insertPoint, lbound, rbound, bestTour, bestCost, forward);
        if(success) while(doAllOpts(bestTour, doNothing, lbound-tours.back().begin()+bestTour.begin(), rbound-tours.back().end()+bestTour.end(), bestTour, bestCost, forward));
        else return false;
        
        if(!forward)
        {
            begin_offset^=(end_offset^=begin_offset^=end_offset);
            std::cout << "swaping offset " << begin_offset << " " << end_offset << std::endl;
        }
        
        tours.back() = std::move(bestTour);
        attended_points.insert(point);
        std::cout << "best result is " << std::endl << tours.back() << std::endl;
        return true;
    }
    
    bool tryInsertRange(TourBound from, TourBound to, TourBound lbound, TourBound rbound)
    {
        std::cout << "try insert range";
        if(rbound-lbound!=0)
            std::cout << " somewhere between " << *lbound << " and " << *(rbound - 1);
        std::cout << std::endl;
        std::vector<PointID> bestTour;
        int bestCost = UNREACHED_WEIGHT;
        bool forward;
        
        auto insertRange = [from, to](std::vector<PointID> &tour){ tour.insert(tour.end(), from, to); };
        auto doNothing = [](std::vector<PointID> &tour){ };
        bool success = doAllOpts(tours.back(), insertRange, lbound, rbound, bestTour, bestCost, forward);
        if(success) while(doAllOpts(bestTour, doNothing, lbound-tours.back().begin()+bestTour.begin(), rbound-tours.back().end()+bestTour.end(), bestTour, bestCost, forward));
        else return false;
        
        if(!forward)
        {
            begin_offset^=(end_offset^=begin_offset^=end_offset);
            std::cout << "swaping offset " << begin_offset << " " << end_offset << std::endl;
        }
        
        tours.back() = std::move(bestTour);
        std::cout << "best result is " << std::endl << tours.back() << std::endl;
        return true;
    }
    
    void startRoutingFromCore(CoreID cur_core)
    {
        begin_offset = 0;
        end_offset = 0;
        PointID from_point = 0;
        PointID start_point = core_start_points[cur_core];
        while(true)
        {
            std::cout << "start routing throw " << cur_core << " core" << std::endl;
            std::cout << "begin_offset is " << begin_offset << std::endl;
            std::cout << "end_offset is " << end_offset << std::endl;
            std::cout << "from_point is " << from_point << std::endl;
            std::cout << "start_point is " << start_point << std::endl;
            std::vector<PointID> tempTour;
            coreCWVisitation(tempTour, cur_core, from_point, start_point);
            std::cout << "CWVisitation is" << std::endl << tempTour << std::endl;
            if(tryInsertRange(tempTour.begin(), tempTour.end(), tours.back().begin() + begin_offset, tours.back().end() - end_offset))
            {
                std::cout << "start inserting inner core points" << std::endl;
                attended_cores.insert(cur_core);
                for(PointID point : core_points[cur_core])
                    if(attended_points.find(point) == attended_points.end())
                        tryInsert(point, tours.back().begin() + begin_offset, tours.back().end() - end_offset);
                    
                std::cout << "core solution is" << std::endl << tours.back() << std::endl;
                std::cout << "start inserting tails" << std::endl;
                    
                auto is_core_point=[&](PointID i) { return core_points[cur_core].find(i) != core_points[cur_core].end(); };
                for(PointID point : core_tails[cur_core])
                    if(attended_points.find(point) == attended_points.end() && isSubSet(attended_cores, tail_cores[point]))
                        //tryInsert(point, tours.back().end() - end_offset, tours.back().begin() + begin_offset);
                    {
                        TourBound min_lbound=tours.back().end() - end_offset;
                        TourBound max_rbound=tours.back().begin() + begin_offset;
                        for(PointID gate : tail_gates[point])
                        {
                            auto gate_iter = std::find(tours.back().begin() + begin_offset, tours.back().end() - end_offset, gate);
                            auto rbound = std::find_if(gate_iter + 1, tours.back().end() - end_offset, is_core_point);
                            if(rbound != tours.back().end() - end_offset) ++rbound;
                            auto lbound = std::find_if(std::reverse_iterator<TourBound>(gate_iter-1), 
                                                    std::reverse_iterator<TourBound>(tours.back().begin() + begin_offset), 
                                                    is_core_point).base();
                            if(min_lbound > lbound) min_lbound = lbound;
                            if(max_rbound < rbound) max_rbound = rbound;
                        }
                        tryInsert(point, min_lbound, max_rbound);
                    }
                std::cout << "tails solution is" << std::endl << tours.back() << std::endl;
            }
            
            std::cout << "start looking for next core" << std::endl;
            CoreID next_core, from_core;
            int min_dist = UNREACHED_WEIGHT;
            PointID min_out_gate, min_in_gate;
            for(CoreID core : attended_cores)
                for(PointID cross_tail : core_tails[core])
                    for(PointID out_gate : tail_gates[cross_tail])
                        if(point_cores[out_gate] == core)
                            for(PointID in_gate : tail_gates[cross_tail])
                                if(attended_cores.find(point_cores[in_gate]) == attended_cores.end())
                                    if(min_dist > time_matrix(out_gate, in_gate))
                                    {
                                        min_dist = time_matrix(out_gate, in_gate);
                                        next_core = point_cores[in_gate];
                                        min_out_gate = out_gate;
                                        min_in_gate = in_gate;
                                        from_core = core;
                                    }
            if(min_dist == UNREACHED_WEIGHT)
                break;
            
            std::cout << "best coice is " << next_core << " conected whis " << from_core << " throw " << min_out_gate << "-" << min_in_gate << std::endl;
            std::cout << "start looking for insert bounds" << std::endl;
            
            auto is_core_point=[&](PointID i) { return core_points[from_core].find(i) != core_points[from_core].end(); };
            TourBound min_lbound=tours.back().end() - end_offset;
            TourBound max_rbound=tours.back().begin() + begin_offset;
            for(PointID cross_tail : core_tails[from_core])
                if(tail_cores[cross_tail].find(next_core) != tail_cores[cross_tail].end())
                    for(PointID out_gate : tail_gates[cross_tail])
                        if(point_cores[out_gate] == from_core)
                        {
                            auto gate_iter = std::find(tours.back().begin() + begin_offset, tours.back().end() - end_offset, out_gate);
                            auto rbound = std::find_if(gate_iter + 1, tours.back().end() - end_offset, is_core_point);
                            if(rbound != tours.back().end() - end_offset) ++rbound;
                            auto lbound = std::find_if(std::reverse_iterator<TourBound>(gate_iter-1), 
                                                    std::reverse_iterator<TourBound>(tours.back().begin() + begin_offset), 
                                                    is_core_point).base();
                            if(min_lbound > lbound) min_lbound = lbound;
                            if(max_rbound < rbound) max_rbound = rbound;
                        }
            std::cout << "best coice are " << *min_lbound << " " << *(max_rbound-1) << std::endl;
            begin_offset = min_lbound - tours.back().begin();
            end_offset = tours.back().end() - max_rbound;
            std::cout << "new offset are " << begin_offset << " " << end_offset << std::endl;
            //tryInsert(min_out_gate, min_lbound, max_rbound);
            cur_core = next_core;
            from_point = min_out_gate;
            start_point = min_in_gate;
        }
    }
    
    template<typename T>
    bool isSubSet(const std::set<T> &set, const std::set<T> &subset)
    {
        typename std::set<T>::const_iterator 
            first1=set.begin(), first2=subset.begin(), 
             last1=set.end(),    last2=subset.end();
        while (first1!=last1 && first2!=last2)
        {
            if (*first1<*first2) ++first1;
            else if (*first2<*first1) return false;
            else { ++first1; ++first2; }
        }
        return first2==last2;
    }
public:
    
    template <class Container, class FixedPointCoordinate>
    GraphLogistic(const unsigned n_points,
                  const ublas::matrix<int> time_matrix,
                  const ublas::matrix<int> length_matrix,
                  const std::vector<Container> full_forward_graph_container,
                  const std::vector<Container> full_reverse_graph_container,
                  const std::vector<FixedPointCoordinate> &coordinates_container)
        : n_points(n_points),
          time_matrix(time_matrix),
          length_matrix(length_matrix),
          nearest_forward_graph(n_points),
          nearest_reverse_graph(n_points),
          core_distances(0)
    {
        int i=0;
        for(const auto &row : full_forward_graph_container) {
            full_forward_graph.emplace_back(row.begin(), row.end());
            if(row.size())
            {
                //int threshold = UNREACHED_WEIGHT;//row[std::min<unsigned>(row.size(), NEAREST_RADIUS) - 1];
                for(const PointID j : row)
                    //if(time_matrix(i, j) <= threshold) 
                    {
                        nearest_forward_graph[i].insert(j);
                        nearest_reverse_graph[j].insert(i);
                        if(time_matrix(i, j)!=0 && nearest_forward_graph[i].size() >= NEAREST_RADIUS)
                            break;
                        //{
                        //    if(full_forward_graph_container[j].size())
                        //        threshold = std::min<int>(1.1 * time_matrix(i, j), time_matrix(i, j) + time_matrix(j, full_forward_graph_container[j].front())/2);
                        //    else threshold = time_matrix(i, j);
                        //}
                    }
            }
            ++i;
        }
        for(const auto &row : full_reverse_graph_container)
            full_reverse_graph.emplace_back(row.begin(), row.end());
        for(auto &coor : coordinates_container)
            coordinates.emplace_back(coor.lat/1000000.0, coor.lon/1000000.0);
    }
    
    void run()
    {
        //findChains();
        findCorePoints();
        //disunitCores();
        findCoreForwardTails();
        findCoreBackwardTails();
        linkUnattachedCores();
        
        tours.emplace_back();
        attended_points.insert(0);
        auto max=std::max_element(core_distances.begin(), core_distances.end());
        //while(!cores_order.empty())
        //{
        if(n_cores)
            startRoutingFromCore(std::distance(core_distances.begin(), max));
            //cores_order.erase(--cores_order.end());
        //}
        
    }
    
    void render(std::vector<char> &output)
    {
        JSON::Object json_root;
        JSON::Array tours_array;
        JSON::Array nearest_array;
        for(auto tour : tours)
        {
            JSON::Array tour_points_array;
            tour_points_array.values.push_back(0);
            for(PointID p : tour)
                tour_points_array.values.push_back(p);
            tour_points_array.values.push_back(0);
            tours_array.values.push_back(tour_points_array);
        }
        json_root.values["tours"] = tours_array;
        for(PointID start=0; start<n_points; ++start)
        {
            JSON::Array nearest_array_row;
            for(const PointID point : nearest_forward_graph[start])
                nearest_array_row.values.push_back(point);
            nearest_array.values.push_back(nearest_array_row);
        }
        json_root.values["nearest"] = nearest_array;
        json_root.values["n"] = n_points;
        JSON::render(output, json_root);
        
    }
};
    
#endif //GRAPG_LOGISTIC_H