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

typedef unsigned PointID;
typedef unsigned CoreID;
typedef unsigned ChainID;
typedef std::pair<double, double> Coordinate;

template<typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
    std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, " "));
    return out;
}

const unsigned SPECIAL_ID = std::numeric_limits<unsigned>::max();

class ChainBase {
    
    ChainID id;
    CoreID coreInner;
    bool is_attend;
    std::set<CoreID> cores_data;
    std::set<ChainID> gates_data;
    std::vector<ChainID> &chains;
    std::vector<CoreID> &cores;
public:
    ChainBase(ChainID id, std::vector<ChainID> &chains, std::vector<CoreID> &cores)
        : id(id), coreInner(SPECIAL_ID), is_attend(false), chains(chains), cores(cores)
    {
    }
    
    void setInner(CoreID core)
    {
        coreInner = core;
        cores_data.insert(core);
        gates_data.insert(id);
    }
    
    bool isUnattached()
    {
        return cores_data.empty();
    }
    
    const std::set<CoreID>& getCores() const
    {
        return cores_data;
    }
    
    void bindTail(ChainID to)
    {
        for(const CoreID &c : chains[to].getCores())
        {
            cores[c].insert(id);
            cores_data.insert(c);
        }
        gates_data.insert(chains[to].gates_data.begin(), chains[to].gates_data.end());
    }
    
    void bindTail(CoreID core, ChainID gate)
    {
        cores[core].insert(id);
        cores_data.insert(core);
        gates_data.insert(gate);
    }
    
    bool isInner() const
    {
        return coreInner != SPECIAL_ID;
    }
    
    bool isInner(CoreID core) const
    {
        return coreInner == core;
    }
    
    bool sameCoreInner(ChainID with) const
    {
        return isInner() && coreInner==chains[with].coreInner;
    }
    
    virtual void attend()
    {
        is_attend = true;
    }
    
    bool isAttend() const
    {
        return is_attend;
    }
    
    virtual ChainID reverse() const
    {
        return id;
    }
    
    virtual Coordinate& getCoordinates() const = 0;
    
    virtual Coordinate& getCoordinatesOut()
    {
        return getCoordinates();
    }
};

class Point : ChainBase {
    PointID point;
};

class OnedirChain : ChainBase {
    std::vector<PointID> points;
};

class ForwardChain : ChainBase {
    friend ReverseChain;
    ChainID reverse;
    std::vector<ChainID> points;
};

class ReverseChain : ChainBase {
    ChainID base;
};

class Core {
    CoreID id;
    std::set<ChainID> inners_data, tails_data;
    std::vector<ChainID> &chains;
    std::vector<CoreID> &cores;
public:
    void insert(ChainID inner)
    {
        inners_data.insert(inner);
    }
    unsigned size() const
    {
        return inners_data.size();
    }
    const std::set<ChainID>& getInners() const
    {
        return inners_data;
    }
    const std::set<ChainID>& getTails() const
    {
        return tails_data;
    }
};
class Tour {
    std::vector<ChainID> data;
    std::map<ChainID, unsigned> indexes;
public:
    typedef typename std::vector<ChainID>::iterator iterator;
    iterator begin(unsigned lbound = 0) const
    {
        return data.begin() + lbound;
    }
    iterator end(unsigned rbound = 0) const
    {
        return data.end() + rbound;
    }
    unsigned size()
    {
        return data.size();
    }
    int cost(){
        return 0;
    }
    void copy(iterator begin, iterator end)
    {
        std::for_each(begin, end, this->add);
    }
    void reverse_copy(iterator begin, iterator end)
    {
        std::for_each(std::reverse_iterator<iterator>(end), std::reverse_iterator<iterator>(begin), [&](ChainID chain){
            this->add(chains[chain].reverse());
        });
    }
    void findBounds(unsigned int begin_offset, unsigned int end_offset, ChainID gate, unsigned int lbound, unsigned int rbound)
    {
        
    }
    void add(ChainID cur_node)
    {
        indexes.insert(chain, data.size());
        data.push_back(chain);
    }
    static Tour make_new(Tour templateTour, int add_count)
    {
        Tour newTour;
        newTour.data.reserve(templateTour.data.size() + add_count);
        return newTour;
    }
};

template<typename T>
class Graph {
    unsigned n;
    std::vector<std::set<T>> forward_data, reverse_data;
    
public:
    Graph() : n(0) { }
    Graph(unsigned n) : n(n), forward_data(n), reverse_data(n) { }
    
    void asign(unsigned n)
    {
        forward_data.asign(n, std::set<T>());
        reverse_data.asign(n, std::set<T>());
    }
    
    void insert(T from, T to)
    {
        forward_data[from].insert(to);
        reverse_data[to].insert(from);
    }
    
    template<typename Iterator>
    void insert(T from, Iterator begin, Iterator end)
    {
        std::for_each(begin, end, [&] (T to) {
            forward_data[from].insert(to);
            reverse_data[to].insert(from);
        });
    }
    
    const std::set<T>& forward(T from) const
    {
        return forward_data[from];
    }
    
    const std::set<T>& reverse(T from) const
    {
        return reverse_data[from];
    }
    
    bool has(T from, T to) const
    {
        return forward_data[from].find(to) != forward_data[from].end();
    }
};


class GraphLogistic
{
    typedef std::vector<ChainID>::iterator TourBound;
    typedef std::shared_ptr<std::pair<TourBound, TourBound>> TourBoundsPair;
    
    const int UNREACHED_WEIGHT = std::numeric_limits<int>::max();
    const unsigned NEAREST_RADIUS = 2;
    const unsigned MIN_CORE_RADIUS = 3;
    const unsigned MOVING_MEAN_RADIUS = 0;
    const float MOVING_MEAN_FACTOR = 1.9;
    const unsigned MAX_TOUR_LENGTH = 20;
    
    unsigned n_points, n_cores, n_chains;
    ublas::matrix<int> time_matrix, length_matrix;
    std::vector<Coordinate> coordinates;
    Graph<PointID> full_point_graph;
    Graph<ChainID> full_graph, nearest_graph;
    int tarjan_time;
    std::list<ChainID> tarjan_stack;
    std::vector<int> tarjan_lowlink;
    std::vector<bool> dfs_visited;
    std::vector<Core> cores;
    std::vector<ChainBase> chains;
    std::list<Tour> tours;
    //std::vector<std::set<ChainID>> core_points, core_tails;
    //std::vector<int> core_distances;
    //std::vector<ChainID> core_start_points;
    //std::map<ChainID, CoreID> point_cores;
    //std::vector<std::set<CoreID>> tail_cores;
    //std::vector<std::set<ChainID>> tail_gates;
    //std::map<ChainID, Chain> chains_index;
    //std::set<ChainID> attended_points;
    //std::set<CoreID> attended_cores;
    unsigned begin_offset, end_offset;
    
    void findCorePoints() 
    {
        tarjan_time = 0;
        tarjan_stack.clear();
        dfs_visited.assign(n_points, false);
        tarjan_lowlink.assign(n_points, 0);
        //chains.assign(n_points, Chain());
        cores.clear();
        cores.reserve(n_points * 0.6);

        for (ChainID u = 0; u < n_points; u++)
            if (!dfs_visited[u])
                tarjanDFS(u);
        
        n_cores = cores.size();
    }
    
    void tarjanDFS(ChainID u) 
    {
        tarjan_lowlink[u] = tarjan_time++;
        dfs_visited[u] = true;
        tarjan_stack.push_back(u);
        bool isComponentRoot = true;

        for (ChainID v : nearest_graph.forward(u)) {
            if (!dfs_visited[v])
                tarjanDFS(v);
            if (tarjan_lowlink[u] > tarjan_lowlink[v]) {
                tarjan_lowlink[u] = tarjan_lowlink[v];
                isComponentRoot = false;
            }
        }

        if (isComponentRoot) {
            Core core;
            //std::set<ChainID> core;
            while (true) {
                ChainID x = tarjan_stack.back();
                tarjan_stack.pop_back();
                if(x != 0) core.insert(x);
                tarjan_lowlink[x] = std::numeric_limits<int>::max();
                if (x == u)
                    break;
            }
            if(core.size() >= MIN_CORE_RADIUS)
            {
                cores.push_back(core);
                const CoreID cur_core = cores.size() - 1;
                //int min_dist = UNREACHED_WEIGHT;
                //ChainID min;
                for(ChainID v : core)
                    chains[v].setInner(cur_core);
                //{
                    //point_cores[p] = cur_core;
                    //tail_cores[p].insert(cur_core);
                    //tail_gates[p].insert(p);
                    //tails_dfs_visited[p] = true;
                    //if(min_dist > time_matrix(0, p))
                    //{
                    //    min_dist = time_matrix(0, p);
                    //    min = p;
                    //}
                //}
                //core_distances.push_back(min_dist);
                //core_start_points.push_back(min);
                //cores_order.insert(cur_core);
            }
        }
    }
    
    void findCoreForwardTails()
    {
        dfs_visited.assign(n_points, false);
        //core_tails.assign(n_cores, std::set<ChainID>());

        for (ChainID u = 0; u < n_points; u++)
            if (!dfs_visited[u])
                forwardTailsDFS(u);
        
        for(ChainID u=0; u<n_chains; u++)
            if(chains[u].isUnattached())
                rebindTail(u);
            
    }

    void forwardTailsDFS(ChainID u) 
    {
        dfs_visited[u] = true;
        for (ChainID v : nearest_graph.forward(u)) {
            if (!dfs_visited[v] && !chains[v].isInner())
                forwardTailsDFS(v);
            if(!chains[u].sameCoreInner(v))
                chains[u].bindTail(v);
            /*for(const CoreID &c : tail_cores[v])
            {
                core_tails[c].insert(u);
                tail_cores[u].insert(c);
            }
            if(point_cores.find(v) == point_cores.end() || point_cores.find(u) == point_cores.end())
                tail_gates[u].insert(tail_gates[v].begin(), tail_gates[v].end());*/
        }
    }
    
    void rebindTail(ChainID u)
    {
        int min_w = UNREACHED_WEIGHT;
        ChainID min = SPECIAL_ID;
        for(const ChainID v : full_graph.forward(u))
            if(!chains[v].isUnattached() && min_w > time_matrix(u, v))
            {
                min_w = time_matrix(u, v);
                min = v;
            }
        if(min != SPECIAL_ID)
            chains[u].bindTail(min);
        /*{
            
            for(const CoreID &c : tail_cores[min])
            {
                core_tails[c].insert(u);
                tail_cores[u].insert(c);
            }
            tail_gates[u].insert(tail_gates[min].begin(), tail_gates[min].end());
        }*/
    }
    
    void findCoreBackwardTails()
    {
        dfs_visited.assign(n_points, false);
        for (CoreID c = 0; c < n_cores; c++)
            for(const ChainID g : cores[c].getInners())
                for(const ChainID v : nearest_graph.forward(g))
                    if(!dfs_visited[v] && !chains[v].isInner(c))
                        backwardTailsDFS(c, g, v);
    }

    void backwardTailsDFS(const CoreID c, const ChainID g, const ChainID u) 
    {
        dfs_visited[u] = true;
        chains[u].bindTail(c, g);
        /*core_tails[c].insert(u);
        tail_cores[u].insert(c);
        tail_gates[u].insert(g);*/
        
        if (!chains[u].isInner())
            for (const ChainID v : nearest_graph.forward(u))
                if (!dfs_visited[v] && !chains[v].isInner())
                    backwardTailsDFS(c, g, v);
    }

    void linkUnattachedCores()
    {
        if(n_cores > 1)
            for(CoreID c = 0; c < n_cores; ++c)
            {
                bool unattached = true;
                for(ChainID t : cores[c].getTails())
                    if(chains[t].getCores().size() > 1) {
                        unattached = false;
                        break;
                    }
                if(unattached)
                {
                    ChainID min_from = SPECIAL_ID, min_to = SPECIAL_ID;
                    int min_weight = UNREACHED_WEIGHT;
                    for(ChainID from : cores[c].getTails())
                        for(ChainID to : full_graph.forward(from))
                            if(!chains[to].isInner(c) && min_weight > time_matrix(from, to))
                            {
                                min_weight = time_matrix(from, to);
                                min_from = from;
                                min_to = to;
                            }
                    chains[min_to].bindTail(min_from);
                    min_weight = UNREACHED_WEIGHT
                    for(ChainID to : cores[c].getTails())
                        for(ChainID from : full_graph.reverse(to))
                            if(!chains[from].isInner(c) && min_weight > time_matrix(from, to))
                            {
                                min_weight = time_matrix(from, to);
                                min_from = from;
                                min_to = to;
                            }
                    chains[min_to].bindTail(min_from);
                    /*tail_cores[min_to].insert(c);
                    tail_gates[min_to].insert(tail_gates[min_from].begin(), tail_gates[min_from].end());
                    tail_gates[min_from].insert(tail_gates[min_to].begin(), tail_gates[min_to].end());*/
                }
            }
    }
    
    float vectprod(Coordinate &A, Coordinate &B, Coordinate &C)
    { 
        return (B.first-A.first) * (C.second-B.second) - (B.second-A.second) * (C.first-B.first); 
    }
    
    void coreCWVisitation(Tour &tour, CoreID cur_core, ChainID from, ChainID first)
    {
        ChainID second;
        ChainID cur_node = first;
        ChainID prev_node = from;
        //std::set<ChainID> attended_points;
        std::set<std::pair<ChainID, ChainID>> attended_edges;
        chains[cur_node].attend();
        //attended_points.insert(cur_node);
        tour.add(cur_node);
        while(true)
        {
            auto i = full_graph.forward(cur_node).begin();
            while(!chains[*i].isInner(cur_core)) ++i;
            ChainID next=*i;
            bool leftmost=vectprod(chains[prev_node].getCoordinatesOut(), chains[cur_node].getCoordinates(), chains[next].getCoordinates())<0;
            for(++i; i!=full_graph.forward(cur_node).end(); ++i)
            {
                if(!chains[*i].isInner(cur_core)) continue;
                if(!leftmost && vectprod(chains[prev_node].getCoordinatesOut(), chains[cur_node].getCoordinates(), chains[*i].getCoordinates())<0){
                    next=*i;
                    leftmost=true;
                }
                else if(*i!=prev_node && vectprod(chains[cur_node].getCoordinates(), chains[next].getCoordinates(), chains[*i].getCoordinates())<0
                    && (!leftmost || vectprod(chains[prev_node].getCoordinatesOut(), chains[cur_node].getCoordinates(), chains[*i].getCoordinates())<0))
                        next=*i;
            }
            prev_node=cur_node;
            cur_node=next;
            auto cur_edge = std::make_pair(prev_node, cur_node);
            if(attended_edges.find(cur_edge) == attended_edges.end())
                attended_edges.insert(cur_edge);
            else break;
            if(!chains[cur_node].isAttend())
            {
                chains[cur_node].attend();
                //attended_points.insert(cur_node);
                tour.add(cur_node);
            }
            //if(tours.back().size() >= MAX_TOUR_LENGTH)
            //    tours.emplace_back();
        }
    }
    
    /*int cost(const std::vector<ChainID> &tour)
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
    }*/
    
    bool checkOpt(const Tour &tempTour, Tour &bestTour)
    {
        //int tempCost = cost(tempTour);
        if(bestTour > tempTour)
        {
            std::cout << "route " << tempTour << " cost " << tempTour.cost() << " is best!" << std::endl;
            bestTour = std::move(tempTour);
            //bestCost = tempCost;
            return true;
        }
        return false;
    }
    
    template<class Action> 
    void do1opt(const Tour &tour, Action action, int action_add_count, unsigned lbound, unsigned rbound, Tour &bestTour, bool &forward)
    {
        std::cout << "opt1" << std::endl; 
        for(auto iter = tour.begin(lbound); iter <= tour.end(rbound); ++iter)
        {
            Tour tempTour = Tour::make_new(tour, action_add_count);
            
            tempTour.copy(tour.begin(), iter);
            action(tempTour);
            tempTour.copy(iter, tour.end());
            
            if(checkOpt(tempTour, bestTour))
                forward = true;
        }
    }
    
    template<class Action> 
    void do2opt(const Tour &tour, Action action, int action_add_count, unsigned lbound, unsigned rbound, Tour &bestTour, bool &forward)
    {
        std::cout << "opt2" << std::endl;
        for(auto iter1=tour.begin(lbound); iter1 <= tour.end(rbound); ++iter1)
            for(auto iter2=tour.begin(lbound); iter2 <= tour.end(rbound); ++iter2)
            {
                Tour tempTour = Tour::make_new(tour, action_add_count);
                
                if(iter2 > iter1 + 1)
                {
                    tempTour.copy(tour.begin(), iter1);
                    action(tempTour);
                    tempTour.reverse_copy(iter1, iter2);
                    tempTour.copy(iter2, tour.end());
                }
                else if(iter2 < iter1)
                {
                    tempTour.reverse_copy(iter1, tour.end());
                    tempTour.copy(iter2, iter1);
                    action(tempTour);
                    tempTour.reverse_copy(tour.begin(), iter2);
                }
                else continue;
                
                if(checkOpt(tempTour, bestTour))
                    forward = iter2 > iter1 + 1;
            }
    }
    
    template<class Action> 
    void do3opt(const Tour &tour, Action action, int action_add_count, unsigned lbound, unsigned rbound, Tour &bestTour, bool &forward)
    {
        std::cout << "opt3" << std::endl;
        for(auto iter1=tour.begin(lbound); iter1 <= tour.end(rbound); ++iter1)
            for(auto iter2=iter1 + 2; iter2 <= tour.end(rbound); ++iter2)
                for(auto iter3=tour.begin(lbound); iter3 <= tour.end(rbound); ++iter3)
                {
                    Tour tempTour = Tour::make_new(tour, action_add_count);
                    
                    if(iter3 > iter2 + 1)
                    {
                        tempTour.copy(tour.begin(), iter1);
                        action(tempTour);
                        tempTour.reverse_copy(iter1, iter2);
                        tempTour.reverse_copy(iter2, iter3);
                        tempTour.copy(iter3, tour.end());
                        
                    }
                    else if(iter3 < iter1)
                    {
                        tempTour.reverse_copy(iter2, tour.end());
                        tempTour.copy(iter3, iter1);
                        action(tempTour);
                        tempTour.reverse_copy(iter1, iter2);
                        tempTour.reverse_copy(tour.begin(), iter3);
                    }
                    else continue;
                    
                    if(checkOpt(tempTour, bestTour))
                        forward = iter3 > iter2 + 1;
                }
    }
    
    template<class Action> 
    bool doAllOpts(const Tour &tour, Action action, int action_add_count, unsigned lbound, unsigned rbound, Tour &bestTour, bool &forward)
    {
        const int prevCost = bestTour.cost();
        do1opt(tour, action, action_add_count, lbound, rbound, bestTour, forward);
        do2opt(tour, action, action_add_count, lbound, rbound, bestTour, forward);
        do3opt(tour, action, action_add_count, lbound, rbound, bestTour, forward);
        return bestTour.cost() < prevCost;
    }

    bool tryInsert(ChainID point, unsigned lbound, unsigned rbound)
    {
        std::cout << "try insert " << point << " somewhere between " << lbound << " and " << rbound << std::endl;
        Tour bestTour;
        bool forward = true;
        
        auto insertPoint = [point](Tour &tour){ tour.add(point); };
        bool success = doAllOpts(tours.back(), insertPoint, 1, lbound, rbound, bestTour, forward);
        if(!success) return false;
        
        if(!forward)
        {
            begin_offset^=(end_offset^=begin_offset^=end_offset);
            std::cout << "swaping offset " << begin_offset << " " << end_offset << std::endl;
        }
        
        
        tours.back() = std::move(bestTour);
        chains[point].attend();
        std::cout << "best result is " << std::endl << tours.back() << std::endl;
        return true;
    }
    
    bool tryInsertRange(TourBound from, TourBound to, unsigned lbound, unsigned rbound)
    {
        std::cout << "try insert range somewhere between " << lbound << " and " << rbound << std::endl;
        Tour bestTour;
        bool forward = true;
        
        auto insertRange = [from, to](Tour &tour){ tour.copy(from, to); };
        bool success = doAllOpts(tours.back(), insertRange, std::distance(from, to), lbound, rbound, bestTour, forward);
        if(!success) return false;
        
        if(!forward)
        {
            begin_offset^=(end_offset^=begin_offset^=end_offset);
            std::cout << "swaping offset " << begin_offset << " " << end_offset << std::endl;
        }
        
        tours.back() = std::move(bestTour);
        std::for_each(from, to, [&] (ChainID point) { chains[point].attend(); });
        std::cout << "best result is " << std::endl << tours.back() << std::endl;
        return true;
    }
    
    bool postoptimize(unsigned lbound, unsigned rbound)
    {
        std::cout << "try postoptimize somewhere between " << lbound << " and " << rbound << std::endl;
        Tour bestTour = tours.back();
        bool forward = true;
        
        auto doNothing = [](Tour &tour){ };
        while(doAllOpts(bestTour, doNothing, 0, lbound, rbound, bestTour, forward)) 
            if(!forward)
            {
                begin_offset^=(end_offset^=begin_offset^=end_offset);
                std::cout << "swaping offset " << begin_offset << " " << end_offset << std::endl;
            }
        if(tours.back() >= bestTour) return false;
        
        tours.back() = std::move(bestTour);
        std::cout << "best result is " << std::endl << tours.back() << std::endl;
        return true;
    }
    
    void startRoutingFromCore(CoreID cur_core)
    {
        begin_offset = 0;
        end_offset = 0;
        ChainID from_point = 0;
        ChainID start_point = core_start_points[cur_core];
        while(true)
        {
            std::cout << "start routing throw " << cur_core << " core" << std::endl;
            std::cout << "begin_offset is " << begin_offset << std::endl;
            std::cout << "end_offset is " << end_offset << std::endl;
            std::cout << "from_point is " << from_point << std::endl;
            std::cout << "start_point is " << start_point << std::endl;
            std::vector<ChainID> tempTour;
            coreCWVisitation(tempTour, cur_core, from_point, start_point);
            std::cout << "CWVisitation is" << std::endl << tempTour << std::endl;
            if(tryInsertRange(tempTour.begin(), tempTour.end(), begin_offset, end_offset))
            {
                std::cout << "start inserting inner core points" << std::endl;
                cores[cur_core].attend();
                for(ChainID point : cores[cur_core].getInners())
                    if(!chains[point].isAttend())
                        tryInsert(point, begin_offset, end_offset);
                    
                std::cout << "core solution is" << std::endl << tours.back() << std::endl;
                std::cout << "start inserting tails" << std::endl;
                    
                //auto is_core_point=[&](ChainID i) { return chains[i].isInner(cur_core); };
                for(ChainID point : cores[cur_core].getTails())
                    if(!chains[point].isAttend() && 
                        std::all_of(chains[point].getCores().begin(), 
                                    chains[point].getCores().end(), 
                                    [&](CoreID c){ return cores[c].isAttend(); }))
                        //tryInsert(point, tours.back().end() - end_offset, tours.back().begin() + begin_offset);
                    {
                        unsigned min_lbound=end_offset;
                        unsigned max_rbound=begin_offset;
                        for(ChainID gate : chains[point].getGates())
                        {
                            unsigned lbound, rbound;
                            tours.back().findBounds(begin_offset, end_offset, gate, lbound, rbound);
                            /*auto gate_iter = std::find(tours.back().begin() + begin_offset, tours.back().end() - end_offset, gate);
                            auto rbound = std::find_if(gate_iter + 1, tours.back().end() - end_offset, is_core_point);
                            if(rbound != tours.back().end() - end_offset) ++rbound;
                            auto lbound = std::find_if(std::reverse_iterator<TourBound>(gate_iter-1), 
                                                    std::reverse_iterator<TourBound>(tours.back().begin() + begin_offset), 
                                                    is_core_point).base();*/
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
            ChainID min_out_gate, min_in_gate;
            for(CoreID core=0; core<n_cores; ++core)
                if(cores[core].isAttend())
                    for(ChainID cross_tail : cores[core].getTails())
                        for(ChainID out_gate : chains[cross_tail].getGates())
                            if(chains[out_gate].isInner(core))
                                for(ChainID in_gate : chains[cross_tail].getGates())
                                    if(!cores[chains[in_gate].getCoreInner()].isAttend())
                                        if(min_dist > time_matrix(out_gate, in_gate))
                                        {
                                            min_dist = time_matrix(out_gate, in_gate);
                                            next_core = chains[in_gate].getCoreInner();
                                            min_out_gate = out_gate;
                                            min_in_gate = in_gate;
                                            from_core = core;
                                        }
            if(min_dist == UNREACHED_WEIGHT)
                break;
            
            std::cout << "best coice is " << next_core << " conected whis " << from_core << " throw " << min_out_gate << "-" << min_in_gate << std::endl;
            std::cout << "start looking for insert bounds" << std::endl;
            

            unsigned min_lbound=end_offset;
            unsigned max_rbound=begin_offset;
            for(ChainID cross_tail : cores[from_core].getTails())
                if(chains[cross_tail].hasCore(next_core))
                    for(ChainID out_gate : chains[cross_tail].getGates())
                        if(chains[out_gate].isInner(from_core))
                        {
                            unsigned lbound, rbound;
                            tours.back().findBounds(begin_offset, end_offset, out_gate, lbound, rbound);
                            if(min_lbound > lbound) min_lbound = lbound;
                            if(max_rbound < rbound) max_rbound = rbound;
                        }
            begin_offset = min_lbound;
            end_offset = max_rbound;
            cur_core = next_core;
            from_point = min_out_gate;
            start_point = min_in_gate;
            std::cout << "new offset are " << begin_offset << " " << end_offset << std::endl;
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
          //nearest_forward_graph(n_points),
          //nearest_reverse_graph(n_points),
          //core_distances(0)
    {
        /*int i=0;
        for(const auto &row : full_forward_graph_container) {
            full_forward_graph.emplace_back(row.begin(), row.end());
            if(row.size())
            {
                //int threshold = UNREACHED_WEIGHT;//row[std::min<unsigned>(row.size(), NEAREST_RADIUS) - 1];
                for(const ChainID j : row)
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
            coordinates.emplace_back(coor.lat/1000000.0, coor.lon/1000000.0);*/
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
        chains[0].attend();
        /*auto max=std::max_element(core_distances.begin(), core_distances.end());
        //while(!cores_order.empty())
        //{
        if(n_cores)
            startRoutingFromCore(std::distance(core_distances.begin(), max));
            //cores_order.erase(--cores_order.end());
        //}*/
        
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
            for(ChainID p : tour)
                tour_points_array.values.push_back(p);
            tour_points_array.values.push_back(0);
            tours_array.values.push_back(tour_points_array);
        }
        json_root.values["tours"] = tours_array;
        for(ChainID start=0; start<n_points; ++start)
        {
            JSON::Array nearest_array_row;
            for(const ChainID point : nearest_graph.forward(start))
                nearest_array_row.values.push_back(point);
            nearest_array.values.push_back(nearest_array_row);
        }
        json_root.values["nearest"] = nearest_array;
        json_root.values["n"] = n_points;
        JSON::render(output, json_root);
        
    }
};
    
#endif //GRAPG_LOGISTIC_H