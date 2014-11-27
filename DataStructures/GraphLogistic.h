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
    
    class CoreByDistanseComparator
    {
        std::vector<int> &core_distances;
        
    public:
        CoreByDistanseComparator(std::vector<int> &core_distance)
            : core_distances(core_distances) { 
                
            }
            
        bool operator()(const CoreID a, const CoreID b) const
        { return core_distances[a] < core_distances[b]; }
    };
    
    const int UNREACHED_WEIGHT = std::numeric_limits<int>::max();
    const unsigned NEAREST_RADIUS = 2;
    const unsigned MIN_CORE_RADIUS = 3;
    const unsigned MOVING_MEAN_RADIUS = 0;
    const float MOVING_MEAN_FACTOR = 1.5;
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
    std::vector<std::set<PointID>> core_points, core_tails, core_gates;
    std::vector<int> core_distances;
    std::vector<PointID> core_start_points;
    std::set<CoreID, CoreByDistanseComparator> cores_order;
    std::vector<std::map<CoreID, std::set<PointID>>> tail_cores;
    std::set<PointID> cored_points;
    std::map<ChainID, Chain> chains_index;
    std::list<std::vector<PointID>> tours;
    std::map<PointID, std::pair<PointID, PointID>> tail_gate_bounds;
    std::set<PointID> attended_points;
    unsigned begin_offset, end_offset;
    
    void findCorePoints() 
    {
        tarjan_time = 0;
        tarjan_stack.clear();
        tarjan_dfs_visited.assign(n_points, false);
        tails_dfs_visited.assign(n_points, false);
        tarjan_lowlink.assign(n_points, 0);
        tail_cores.assign(n_points, std::map<CoreID, std::set<PointID>>());
        core_points.clear();
        core_points.reserve(n_points * 0.6);

        for (PointID u = 0; u < n_points; u++)
            if (!tarjan_dfs_visited[u])
                tarjan_dfs(u);
        
        n_cores = core_points.size();
    }
    
    void findCoreTails()
    {
        core_tails.assign(n_cores, std::set<PointID>());

        for (PointID u = 0; u < n_points; u++)
            if (!tails_dfs_visited[u])
                tails_dfs(u);
        
        auto &tail_cores_ref=tail_cores;
        rebindTail(boost::make_counting_iterator(0u), 
                   boost::make_counting_iterator(n_points),
                   [&tail_cores_ref](PointID p) { return tail_cores_ref[p].empty(); });
            
    }
    
    template<class Iter, class Filter>
    void rebindTail(Iter begin, Iter end, Filter filter)
    {
        for(Iter p=begin; p<end; p++)
            if(filter(*p))
            {
                for(const auto &c : tail_cores[*p])
                    core_tails[c.first].erase(*p);
                tail_cores[*p].clear();
                int min_w = UNREACHED_WEIGHT;
                PointID min;
                for(const PointID i : full_forward_graph[*p])
                    if(!filter(i) && min_w > time_matrix(*p, i))
                    {
                        min_w = time_matrix(*p, i);
                        min = i;
                    }
                if(min_w != UNREACHED_WEIGHT)
                    for(const auto &c : tail_cores[min])
                    {
                        core_tails[c.first].insert(*p);
                        tail_cores[*p][c.first].insert(c.second.begin(), c.second.end());
                    }
            }
    }
    
    void findCoreGates()
    {
        core_gates.assign(n_cores, std::set<PointID>());
        for (CoreID c = 0; c < n_cores; c++)
            for(const PointID u : core_points[c])
                for(const PointID v : nearest_forward_graph[u])
                    if(core_points[c].find(v) == core_points[c].end())
                        core_gates[c].insert(v);
    }

    void tarjan_dfs(PointID u) 
    {
        tarjan_lowlink[u] = tarjan_time++;
        tarjan_dfs_visited[u] = true;
        tarjan_stack.push_back(u);
        bool isComponentRoot = true;

        for (PointID v : nearest_forward_graph[u]) {
            if (!tarjan_dfs_visited[v])
                tarjan_dfs(v);
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
                core.insert(x);
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
                    cored_points.insert(p);
                    tail_cores[p][cur_core].insert(p);
                    tails_dfs_visited[p] = true;
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
      
    void tails_dfs(PointID u) 
    {
        tails_dfs_visited[u] = true;
        for (PointID v : nearest_forward_graph[u]) {
            if (!tails_dfs_visited[v])
                tails_dfs(v);
            for(const auto &c : tail_cores[v])
            {
                core_tails[c.first].insert(u);
                tail_cores[u][c.first].insert(c.second.begin(), c.second.end());
            }
        }
    }

    float vectprod(std::pair<double,double> &A, std::pair<double,double> &B, std::pair<double,double> &C)
    { 
        return (B.first-A.first) * (C.second-B.second) - (B.second-A.second) * (C.first-B.first); 
    }
    
    void coreCWVisitation(std::vector<PointID> &tour, CoreID cur_core, PointID from, PointID first)
    {
        PointID cur_node = first;
        PointID prev_node = from;
        attended_points.insert(cur_node);
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
            if(attended_points.find(next) != attended_points.end())            
                break;
            prev_node=cur_node;
            cur_node=next;
            attended_points.insert(cur_node);
            tour.push_back(cur_node);
            //if(tours.back().size() >= MAX_TOUR_LENGTH)
            //    tours.emplace_back();
        }
    }
    
    int cost(std::vector<PointID> &tour)
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
    
    bool tryInsert(PointID point, TourBound lbound, TourBound rbound)
    {
        std::cout << "try insert " << point << " somewhere between " << *lbound << " and " << *(rbound - 1) << std::endl;
        std::vector<PointID> bestTour;
        int bestCost = UNREACHED_WEIGHT;
        bool forward;
        
        auto insertPoint = [point](std::vector<PointID> &tour){ tour.push_back(point); };
        do1opt(insertPoint, lbound, rbound, bestTour, bestCost, forward);
        do2opt(insertPoint, lbound, rbound, bestTour, bestCost, forward);
        do3opt(insertPoint, lbound, rbound, bestTour, bestCost, forward);
        
        if(bestCost == UNREACHED_WEIGHT) return false;
        
        if(!forward)
            begin_offset^=(end_offset^=begin_offset^=end_offset);
        
        tours.back() = std::move(bestTour);
        std::cout << "best result is " << std::endl << tours.back() << std::endl;
        return true;
    }
    
    bool tryInsertRange(TourBound from, TourBound to, TourBound lbound, TourBound rbound)
    {
        if(rbound-lbound!=0)
            std::cout << "try insert range somewhere between " << *lbound << " and " << *(rbound - 1) << std::endl;
        std::vector<PointID> bestTour;
        int bestCost = UNREACHED_WEIGHT;
        bool forward;
        
        auto insertRange = [from, to](std::vector<PointID> &tour){ tour.insert(tour.end(), from, to); };
        do1opt(insertRange, lbound, rbound, bestTour, bestCost, forward);
        do2opt(insertRange, lbound, rbound, bestTour, bestCost, forward);
        do3opt(insertRange, lbound, rbound, bestTour, bestCost, forward);
        
        if(bestCost == UNREACHED_WEIGHT) return false;
        
        if(!forward)
            begin_offset^=(end_offset^=begin_offset^=end_offset);
        
        tours.back() = std::move(bestTour);
        std::cout << "best result is " << std::endl << tours.back() << std::endl;
        return true;
    }
    
    bool postOptimize(TourBound lbound, TourBound rbound)
    {
        std::cout << "try postoptimize somewhere between " << *lbound << " and " << *(rbound - 1) << std::endl;
        std::vector<PointID> bestTour;
        int bestCost = UNREACHED_WEIGHT;
        bool forward;
        
        auto doNothing = [](std::vector<PointID> &tour){ };
        do1opt(doNothing, lbound, rbound, bestTour, bestCost, forward);
        do2opt(doNothing, lbound, rbound, bestTour, bestCost, forward);
        do3opt(doNothing, lbound, rbound, bestTour, bestCost, forward);
        
        if(bestCost == UNREACHED_WEIGHT) return false;
        
        if(!forward)
            begin_offset^=(end_offset^=begin_offset^=end_offset);
        
        tours.back() = std::move(bestTour);
        std::cout << "best result is " << std::endl << tours.back() << std::endl;
        return true;
    }
    
    template<class Action> 
    void do1opt(Action action, TourBound lbound, TourBound rbound, std::vector<PointID> &bestTour, int &bestCost, bool &forward)
    {
        std::cout << "opt1" << std::endl;
        const auto begin = tours.back().begin();
        const auto end = tours.back().end();        
        for(auto iter=lbound; iter <= rbound; ++iter)
        {
            //if(!checkPointsCompilable(tempTour)) continue;
            
            std::vector<PointID> tempTour;
            tempTour.reserve(tours.back().size() + 1);
            
            std::copy(begin, iter, std::back_inserter(tempTour));
            action(tempTour);
            std::copy(iter, end, std::back_inserter(tempTour));
            
            //if(!checkRouteExists(tempTour)) continue;
            int tempCost = cost(tempTour);
            if(bestCost > tempCost)
            {
                std::cout << "route " << tempTour << " cost " << tempCost << " is best!" << std::endl;
                bestTour = std::move(tempTour);
                bestCost = tempCost;
                forward = true;
            }
            //else std::cout << "route " << tempTour << " cost " << tempCost << std::endl;
        }
    }
    
    template<class Action> 
    void do2opt(Action action, TourBound lbound, TourBound rbound, std::vector<PointID> &bestTour, int &bestCost, bool &forward)
    {
        std::cout << "opt2" << std::endl;
        const auto begin = tours.back().begin();
        const auto end = tours.back().end();        
        for(auto iter1=lbound; iter1 <= rbound; ++iter1)
            for(auto iter2=lbound; iter2 <= rbound; ++iter2)
            {
                //if(!checkPointsCompilable(tempTour)) continue;
                
                std::vector<PointID> tempTour;
                tempTour.reserve(tours.back().size() + 1);
                if(iter2 > iter1 + 1)
                {
                    std::copy(begin, iter1, std::back_inserter(tempTour));
                    action(tempTour);
                    std::reverse_copy(iter1, iter2, std::back_inserter(tempTour));
                    std::copy(iter2, end, std::back_inserter(tempTour));
                }
                else if(iter2 < iter1)
                {
                    std::reverse_copy(iter1, end, std::back_inserter(tempTour));
                    std::copy(iter2, iter1, std::back_inserter(tempTour));
                    action(tempTour);
                    std::reverse_copy(begin, iter2, std::back_inserter(tempTour));
                }
                else continue;
                
                //if(!checkRouteExists(tempTour)) continue;
                int tempCost = cost(tempTour);
                if(bestCost > tempCost)
                {
                    std::cout << "route " << tempTour << " cost " << tempCost << " is best!" << std::endl;
                    bestTour = std::move(tempTour);
                    bestCost = tempCost;
                    forward = iter2 > iter1 + 1;
                }
                //else std::cout << "route " << tempTour << " cost " << tempCost << std::endl;
            }
    }
    
    template<class Action> 
    void do3opt(Action action, TourBound lbound, TourBound rbound, std::vector<PointID> &bestTour, int &bestCost, bool &forward)
    {
        std::cout << "opt3" << std::endl;
        const auto begin = tours.back().begin();
        const auto end = tours.back().end();        
        for(auto iter1=lbound; iter1 <= rbound; ++iter1)
            for(auto iter2=iter1 + 2; iter2 <= rbound; ++iter2)
                for(auto iter3=lbound; iter3 <= rbound; ++iter3)
                {
                    //if(!checkPointsCompilable(tempTour)) continue;
                    
                    std::vector<PointID> tempTour;
                    tempTour.reserve(tours.back().size() + 1);
                    if(iter3 > iter2 + 1)
                    {
                        std::copy(begin, iter1, std::back_inserter(tempTour));
                        action(tempTour);
                        std::reverse_copy(iter1, iter2, std::back_inserter(tempTour));
                        std::reverse_copy(iter2, iter3, std::back_inserter(tempTour));
                        std::copy(iter3, end, std::back_inserter(tempTour));
                        
                    }
                    else if(iter3 < iter1)
                    {
                        std::reverse_copy(iter2, end, std::back_inserter(tempTour));
                        std::copy(iter3, iter1, std::back_inserter(tempTour));
                        action(tempTour);
                        std::reverse_copy(iter1, iter2, std::back_inserter(tempTour));
                        std::reverse_copy(begin, iter3, std::back_inserter(tempTour));
                    }
                    else continue;
                    
                    //if(!checkRouteExists(tempTour)) continue;
                    int tempCost = cost(tempTour);
                    if(bestCost > tempCost)
                    {
                        std::cout << "route " << tempTour << " cost " << tempCost << " is best!" << std::endl;
                        bestTour = std::move(tempTour);
                        bestCost = tempCost;
                        forward = iter3 > iter2 + 1;
                    }
                    //else std::cout << "route " << tempTour << " cost " << tempCost << std::endl;
                }
    }

    //void getTailGateBounds()
    void startRoutingFromCore(CoreID cur_core)
    {
        begin_offset = 0;
        end_offset = 0;
        PointID from_point = 0;
        PointID start_point = core_start_points[cur_core];
        while(true)
        {
            std::vector<PointID> tempTour;
            coreCWVisitation(tempTour, cur_core, from_point, start_point);
            std::cout << "initial solution is" << std::endl << tempTour << std::endl;
            tryInsertRange(tempTour.begin(), tempTour.end(), tours.back().begin() + begin_offset, tours.back().end() - end_offset);
            for(PointID point : core_points[cur_core])
                if(attended_points.find(point) == attended_points.end())
                    tryInsert(point, tours.back().begin() + begin_offset, tours.back().end() - end_offset);
            postOptimize(tours.back().begin() + begin_offset, tours.back().end() - end_offset);
                
            std::cout << "core solution is" << std::endl << tours.back() << std::endl;
                
            auto is_core_point=[&](PointID i) { return core_points[cur_core].find(i) != core_points[cur_core].end(); };
            for(PointID point : core_tails[cur_core])
                if(attended_points.find(point) == attended_points.end() && tail_cores[point].size() == 1)
                {
                    TourBound min_lbound=tours.back().end() - end_offset;
                    TourBound max_rbound=tours.back().begin() + begin_offset;
                    for(PointID gate : tail_cores[point][cur_core])
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
            postOptimize(tours.back().begin() + begin_offset, tours.back().end() - end_offset);
            
            CoreID next_core;
            int min_dist = UNREACHED_WEIGHT;
            PointID min_out_gate, min_in_gate;
            for(PointID out_gate : core_gates[cur_core])
                for(auto core : tail_cores[out_gate])
                    for(PointID in_gate : core.second)
                        if(min_dist > time_matrix(out_gate, in_gate))
                        {
                            min_dist = time_matrix(out_gate, in_gate);
                            next_core = core.first;
                            min_out_gate = out_gate;
                            min_in_gate = in_gate;
                        }
            if(min_dist == UNREACHED_WEIGHT)
                break;
            
            TourBound min_lbound=tours.back().end() - end_offset;
            TourBound max_rbound=tours.back().begin() + begin_offset;
            for(PointID out_gate : core_gates[cur_core])
                if(tail_cores[out_gate].find(next_core) != tail_cores[out_gate].end())
                    for(PointID inner_prev : nearest_reverse_graph[out_gate])
                        if(core_points[cur_core].find(inner_prev) != core_points[cur_core].end())
                        {
                            auto gate_iter = std::find(tours.back().begin() + begin_offset, tours.back().end() - end_offset, inner_prev);
                            auto rbound = std::find_if(gate_iter + 1, tours.back().end() - end_offset, is_core_point);
                            if(rbound != tours.back().end() - end_offset) ++rbound;
                            auto lbound = std::find_if(std::reverse_iterator<TourBound>(gate_iter-1), 
                                                    std::reverse_iterator<TourBound>(tours.back().begin() + begin_offset), 
                                                    is_core_point).base();
                            if(min_lbound > lbound) min_lbound = lbound;
                            if(max_rbound < rbound) max_rbound = rbound;
                        }
            begin_offset += min_lbound - tours.back().begin();
            end_offset += tours.back().end() - max_rbound;
            //tryInsert(min_out_gate, min_lbound, max_rbound);
            cur_core = next_core;
            from_point = min_out_gate;
            start_point = min_in_gate;
        }
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
          core_distances(0),
          cores_order(core_distances)
    {
        int i=0;
        for(const auto &row : full_forward_graph_container) {
            full_forward_graph.emplace_back(row.begin(), row.end());
            if(row.size())
            {
                int threshold = UNREACHED_WEIGHT;//row[std::min<unsigned>(row.size(), NEAREST_RADIUS) - 1];
                for(const PointID j : row)
                    if(time_matrix(i, j) <= threshold) 
                    {
                        nearest_forward_graph[i].insert(j);
                        nearest_reverse_graph[j].insert(i);
                        if(nearest_forward_graph[i].size() == NEAREST_RADIUS) {
                            if(full_forward_graph_container[j].size())
                                threshold = std::min<int>(1.1 * time_matrix(i, j), time_matrix(i, j) + time_matrix(j, full_forward_graph_container[j].front()));
                            else threshold = time_matrix(i, j);
                        }
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
        findCoreTails();
        findCoreGates();
        
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