#ifndef GRAPG_LOGISTIC_H
#define GRAPG_LOGISTIC_H

#include "../DataStructures/JSONContainer.h"
#include "Core.h"
#include "Chain.h"
#include "Graph.h"
#include "Tour.h"


class GraphLogistic
{

    
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
    std::vector<std::shared_ptr<ChainBase>> chains;
    std::list<Tour> tours;
    unsigned begin_offset, end_offset;
    
    void add_into_prechain_list(bool back, 
                                std::shared_ptr<std::list<PointID>> &list, 
                                PointID point,
                                std::map<PointID, std::shared_ptr<std::list<PointID>>> &map)
    {
        if(back) list->push_back(point);
        else list->push_front(point);
        map[point] = list;
    }
    
    void add_glue_into_prechain_list(bool back, 
                                     std::shared_ptr<std::list<PointID>> &list, 
                                     PointID point, 
                                     const std::map<PointID, std::shared_ptr<std::list<PointID>>> &glueM, 
                                     std::set<std::shared_ptr<std::list<PointID>>> &glueL,
                                     std::map<PointID, std::shared_ptr<std::list<PointID>>> &map)
    {
        const auto glue_iter = glueM.find(point);
        if(glue_iter != glueM.end())
        {
            for(PointID p : *glue_iter->second)
                add_into_prechain_list(back, list, p, map);
            glueL.erase(glue_iter->second);
        }
        else add_into_prechain_list(back, list, point, map);
    }
    
    void findChains()
    {
        std::set<std::shared_ptr<std::list<PointID>>> glueL, onedirL, bidirL;
        std::map<PointID, std::shared_ptr<std::list<PointID>>> glueM, onedirM, bidirM;
        for(PointID i=1; i<n_points; i++) //TODO depot
            for(PointID j : full_point_graph.forward(i))
                if(time_matrix(i, j) == 0)
                {
                    auto i_map_iter = glueM.find(i);
                    auto j_map_iter = glueM.find(j);
                    if(i_map_iter!=glueM.end() && j_map_iter!=glueM.end())
                    {
                        if(i_map_iter->second == j_map_iter->second) continue;
                        i_map_iter->second->insert(i_map_iter->second->end(), 
                                                   j_map_iter->second->begin(), 
                                                   j_map_iter->second->end());
                        for(PointID k : *j_map_iter->second)
                            glueM[k] = i_map_iter->second;
                        glueL.erase(j_map_iter->second);
                    }
                    else if(i_map_iter!=glueM.end())
                        i_map_iter->second->push_back(j);
                    else if(j_map_iter!=glueM.end())
                        j_map_iter->second->push_front(i);
                    else
                    {
                        auto list = std::make_shared<std::list<PointID>>();
                        list->push_back(i);
                        list->push_back(j);
                        glueL.insert(list);
                        glueM[i]=list;
                        glueM[j]=list;
                    }
                }
        for(auto gl : glueL)
        {
            for(PointID i : *gl)
                full_point_graph.erase(gl->back(), i);
            std::for_each(gl->begin(), --gl->end(), [&](PointID i){
                for(PointID j : full_point_graph.forward(gl->back()))
                    full_point_graph.erase(i, j);
            });
        }
        for(PointID i=1; i<n_points; i++) //TODO depot
            if(full_point_graph.forward(i).size() == 1
                && full_point_graph.reverse(full_point_graph.forward(i).front()).size() == 1)
            {
                PointID j = full_point_graph.forward(i).front();
                auto i_map_iter = onedirM.find(i);
                auto j_map_iter = onedirM.find(j);
                if(i_map_iter!=onedirM.end() && j_map_iter!=onedirM.end())
                {
                    if(i_map_iter->second == j_map_iter->second) continue;
                    for(PointID p : *j_map_iter->second)
                    {
                        add_into_prechain_list(true, i_map_iter->second, p, onedirM);
                        onedirM[p] = i_map_iter->second;
                    }
                    onedirL.erase(j_map_iter->second);
                }
                else if(i_map_iter!=onedirM.end())
                    add_glue_into_prechain_list(true, i_map_iter->second, j, glueM, glueL, onedirM);
                else if(j_map_iter!=onedirM.end())
                    add_glue_into_prechain_list(false, j_map_iter->second, i, glueM, glueL, onedirM);
                else
                {
                    auto list = std::make_shared<std::list<PointID>>();
                    add_glue_into_prechain_list(true, list, i, glueM, glueL, onedirM);
                    add_glue_into_prechain_list(true, list, j, glueM, glueL, onedirM);
                    onedirL.insert(list);
                }
            }
        for(PointID j=1; j<n_points; j++) //TODO depot
            if(full_point_graph.forward(j).size() == 2 && full_point_graph.reverse(j).size() == 2
                && *full_point_graph.forward_set(j).begin() == *full_point_graph.reverse(j).begin()
                && *(--full_point_graph.forward_set(j).end()) == *(--full_point_graph.reverse(j).end()))
            {
                PointID i = full_point_graph.forward(j).front();
                PointID k = full_point_graph.forward(j).back();
                auto i_map_iter = bidirM.find(i);
                auto j_map_iter = bidirM.find(j);
                auto k_map_iter = bidirM.find(k);
                if(i_map_iter!=bidirM.end() && k_map_iter!=bidirM.end())
                {
                    if(i_map_iter->second == k_map_iter->second) continue;
                    bool back = true;
                    if(j_map_iter == bidirM.end())
                    {
                        back = i_map_iter->second->back() == i;
                        add_into_prechain_list(back, i_map_iter->second, j, bidirM);
                    }
                    else back = i_map_iter->second->back() == j;
                    for(PointID p : *k_map_iter->second)
                    {
                        add_into_prechain_list(back, i_map_iter->second, p, bidirM);
                        bidirM[p] = i_map_iter->second;
                    }
                    bidirL.erase(k_map_iter->second);
                }
                else if(i_map_iter!=bidirM.end())
                {
                    bool back = true;
                    if(j_map_iter == bidirM.end())
                    {
                        back = i_map_iter->second->back() == i;
                        add_into_prechain_list(back, i_map_iter->second, j, bidirM);
                    }
                    else back = i_map_iter->second->back() == j;
                    add_glue_into_prechain_list(back, i_map_iter->second, k, glueM, glueL, bidirM);
                }
                else if(k_map_iter!=bidirM.end())
                {
                    bool back = false;
                    if(j_map_iter == bidirM.end())
                    {
                        back = k_map_iter->second->back() == k;
                        add_into_prechain_list(back, k_map_iter->second, j, bidirM);
                    }
                    else back = k_map_iter->second->back() == j;
                    add_glue_into_prechain_list(back, k_map_iter->second, i, glueM, glueL, bidirM);
                }
                else
                {
                    auto list = std::make_shared<std::list<PointID>>();
                    add_glue_into_prechain_list(true, list, i, glueM, glueL, bidirM);
                    add_into_prechain_list(true, list, j, bidirM);
                    add_glue_into_prechain_list(true, list, k, glueM, glueL, bidirM);
                    bidirL.insert(list);
                }
            }
        for(auto gl : glueL)
            onedirL.insert(gl);
        std::set<PointID> alone_points;
        std::map<PointID, ChainID> points_chain_map;
        chains.emplace_back(new Point(0, 0, chains, cores, time_matrix, length_matrix)); //TODO depot
        points_chain_map[0] = 0;
        for(PointID p = 1; p < n_points; ++p) //TODO depot
            alone_points.insert(alone_points.end(), p);
        for(auto list : onedirL)
        {
            ChainID id = chains.size();
            std::shared_ptr<ChainBase> chain(new OnedirChain(id, *list, chains, cores, time_matrix, length_matrix));
            for(PointID p : *list)
            {
                alone_points.erase(p);
                points_chain_map[p] = id;
            }
            points_chain_map[chain->getInPointID()] = id;
            chains.push_back(chain);
            std::cout<<"chain "<<id<<" is: "<<*list<<std::endl;
        }
        for(auto list : bidirL)
        {
            ChainID id_f = chains.size();
            auto ptr = new ForwardChain(id_f, *list, chains, cores, time_matrix, length_matrix);
            std::shared_ptr<ChainBase> chain_f(ptr);
            chains.push_back(chain_f);
            ChainID id_r = chains.size();
            std::shared_ptr<ChainBase> chain_r(new ReverseChain(id_r, *ptr));
            chains.push_back(chain_r);
            for(PointID p : *list)
            {
                alone_points.erase(p);
                points_chain_map[p] = SPECIAL_ID;
            }
            points_chain_map[chain_f->getInPointID()] = id_f;
            points_chain_map[chain_r->getInPointID()] = id_r;
            std::cout<<"chain "<<id_f<<" is: "<<*list<<std::endl;
            std::cout<<"chain "<<id_r<<" is: ";
            std::reverse_copy(list->begin(), list->end(), std::ostream_iterator<PointID>(std::cout, " "));
            std::cout<<std::endl;
        }
        for(PointID p : alone_points)
        {
            ChainID id = chains.size();
            std::shared_ptr<ChainBase> chain(new Point(id, p, chains, cores, time_matrix, length_matrix));
            points_chain_map[p] = id;
            chains.push_back(chain);
            std::cout<<"chain "<<id<<" is: "<<p<<std::endl;
        }
        n_chains = chains.size();
        full_graph.assign(n_chains);
        for(ChainID u = 0; u < n_chains; ++u)
            for(PointID p : full_point_graph.forward(chains[u]->getOutPointID()))
            {
                const ChainID v = points_chain_map[p];
                if(v == SPECIAL_ID) continue;
                full_graph.insert(u, v);
            }
    }
    
    void buildNearestGraph()
    {
        nearest_graph.assign(n_chains);
        for(ChainID i=0; i<n_chains; i++) {
            int threshold = UNREACHED_WEIGHT;//row[std::min<unsigned>(row.size(), NEAREST_RADIUS) - 1];
            for(const ChainID j : full_graph.forward(i))
                if(chains[i]->timeTo(j) <= threshold) 
                {
                    nearest_graph.insert(i, j);
                    if(nearest_graph.forward(i).size() >= NEAREST_RADIUS)
                    {
                        if(!full_graph.forward(j).empty())
                            threshold = std::min<int>(1.1 * chains[i]->timeTo(j), chains[i]->timeTo(j) + chains[j]->timeTo(full_graph.forward(j).front())/2);
                        else threshold = chains[i]->timeTo(j);
                    }
                }
        }
    }
    
    void findCorePoints() 
    {
        tarjan_time = 0;
        tarjan_stack.clear();
        dfs_visited.assign(n_chains, false);
        tarjan_lowlink.assign(n_chains, 0);
        //chains.assign(n_chains, Chain());
        cores.clear();
        cores.reserve(n_chains * 0.6);

        for (ChainID u = 1; u < n_chains; u++) //TODO depot
            if (!dfs_visited[u])
                tarjanDFS(u);
        
        if(cores.empty())
        {
            Core core(chains, cores);
            for(ChainID i=1; i<n_chains; ++i) //TODO depots
            {
                core.insert(i);
                chains[i]->setInner(0);
            }
            cores.push_back(core);
        }
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
            Core core(chains, cores);
            //std::set<ChainID> core;
            while (true) {
                ChainID x = tarjan_stack.back();
                tarjan_stack.pop_back();
                core.insert(x); 
                tarjan_lowlink[x] = std::numeric_limits<int>::max();
                if (x == u)
                    break;
            }
            if(core.size() >= MIN_CORE_RADIUS)
            {
                cores.push_back(core);
                const CoreID cur_core = cores.size() - 1;
                for(ChainID v : core.getInners())
                    chains[v]->setInner(cur_core);
            }
        }
    }
    
    void findCoreForwardTails()
    {
        dfs_visited.assign(n_chains, false);

        for (ChainID u = 0; u < n_chains; u++)
            if (!dfs_visited[u])
                forwardTailsDFS(u);
        
        for(ChainID u=0; u<n_chains; u++)
            if(chains[u]->isUnattached())
                rebindTail(u);
            
    }

    void forwardTailsDFS(ChainID u) 
    {
        dfs_visited[u] = true;
        for (ChainID v : nearest_graph.forward(u)) {
            if (!dfs_visited[v] && !chains[v]->isInner())
                forwardTailsDFS(v);
            if(!chains[u]->sameCoreInner(v))
                chains[u]->bindTail(v);
        }
    }
    
    void rebindTail(ChainID u)
    {
        int min_w = UNREACHED_WEIGHT;
        ChainID min = SPECIAL_ID;
        for(const ChainID v : full_graph.forward(u))
            if(!chains[v]->isUnattached() && min_w > chains[u]->timeTo(v))
            {
                min_w = chains[u]->timeTo(v);
                min = v;
            }
        if(min != SPECIAL_ID)
            chains[u]->bindTail(min);
    }
    
    void findCoreBackwardTails()
    {
        dfs_visited.assign(n_chains, false);
        for (CoreID c = 0; c < n_cores; c++)
            for(const ChainID g : cores[c].getInners())
                for(const ChainID v : nearest_graph.forward(g))
                    if(!dfs_visited[v] && !chains[v]->isInner(c))
                        backwardTailsDFS(c, g, v);
    }

    void backwardTailsDFS(const CoreID c, const ChainID g, const ChainID u) 
    {
        dfs_visited[u] = true;
        chains[u]->bindTail(c, g);
        
        if (!chains[u]->isInner())
            for (const ChainID v : nearest_graph.forward(u))
                if (!dfs_visited[v] && !chains[v]->isInner())
                    backwardTailsDFS(c, g, v);
    }

    void linkUnattachedCores()
    {
        if(n_cores > 1)
            for(CoreID c = 0; c < n_cores; ++c)
            {
                bool unattached = true;
                for(ChainID t : cores[c].getTails())
                    if(chains[t]->getCores().size() > 1) {
                        unattached = false;
                        break;
                    }
                if(unattached)
                {
                    ChainID min_from = SPECIAL_ID, min_to = SPECIAL_ID;
                    int min_weight = UNREACHED_WEIGHT;
                    for(ChainID from : cores[c].getTails())
                        for(ChainID to : full_graph.forward(from))
                            if(!chains[to]->isInner(c) && min_weight > chains[from]->timeTo(to))
                            {
                                min_weight = chains[from]->timeTo(to);
                                min_from = from;
                                min_to = to;
                            }
                    chains[min_to]->bindTail(min_from);
                    min_weight = UNREACHED_WEIGHT;
                    for(ChainID to : cores[c].getTails())
                        for(ChainID from : full_graph.reverse(to))
                            if(!chains[from]->isInner(c) && min_weight > chains[from]->timeTo(to))
                            {
                                min_weight = chains[from]->timeTo(to);
                                min_from = from;
                                min_to = to;
                            }
                    chains[min_to]->bindTail(min_from);
                }
            }
    }
    
    float vectprod(ChainID a, ChainID b, ChainID c)
    { 
        const auto &A = coordinates[chains[a]->getOutPointID()];
        const auto &B = coordinates[chains[b]->getInPointID()];
        const auto &C = coordinates[chains[c]->getInPointID()];
        return (B.first-A.first) * (C.second-B.second) - (B.second-A.second) * (C.first-B.first); 
    }
    
    void coreCWVisitation(Tour &tour, CoreID cur_core, ChainID from, ChainID first)
    {
        //ChainID second;
        ChainID cur_node = first;
        ChainID prev_node = from;
        std::set<std::pair<ChainID, ChainID>> attended_edges;
        chains[cur_node]->attend();
        tour.add(cur_node);
        while(true)
        {
            auto i = full_graph.forward(cur_node).begin();
            while(i!=full_graph.forward(cur_node).end() && !chains[*i]->isInner(cur_core)) ++i;
            if(i==full_graph.forward(cur_node).end()) break;
            ChainID next=*i;
            bool leftmost=vectprod(prev_node, cur_node, next)<0;
            for(++i; i!=full_graph.forward(cur_node).end(); ++i)
            {
                if(!chains[*i]->isInner(cur_core)) continue;
                if(!leftmost && vectprod(prev_node, cur_node, *i)<0){
                    next=*i;
                    leftmost=true;
                }
                else if(*i!=prev_node && vectprod(cur_node, next, *i)<0
                    && (!leftmost || vectprod(prev_node, cur_node, *i)<0))
                        next=*i;
            }
            prev_node=cur_node;
            cur_node=next;
            auto cur_edge = std::make_pair(prev_node, cur_node);
            if(attended_edges.find(cur_edge) == attended_edges.end())
                attended_edges.insert(cur_edge);
            else break;
            if(!chains[cur_node]->isAttend())
            {
                chains[cur_node]->attend();
                tour.add(cur_node);
            }
        }
    }
    
    bool checkOpt(Tour &tempTour, Tour &bestTour)
    {
        //int tempCost = cost(tempTour);
        if(bestTour > tempTour)
        {
            std::cout << "route " << tempTour << " cost " << tempTour.cost() << " is best!" << std::endl;
            bestTour = std::move(tempTour);
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
            Tour tempTour(tour, action_add_count);
            
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
                Tour tempTour(tour, action_add_count);
                
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
                    Tour tempTour(tour, action_add_count);
                    
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
        Tour bestTour(chains, time_matrix);
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
        chains[point]->attend();
        std::cout << "best result is " << std::endl << tours.back() << std::endl;
        return true;
    }
    
    bool tryInsertRange(TourBound from, TourBound to, unsigned lbound, unsigned rbound)
    {
        std::cout << "try insert range somewhere between " << lbound << " and " << rbound << std::endl;
        Tour bestTour(chains, time_matrix);
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
        std::for_each(from, to, [&] (ChainID point) { chains[point]->attend(); });
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
    
    void startRoutingFromCore(CoreID cur_core, ChainID start_point)
    {
        begin_offset = 0;
        end_offset = 0;
        ChainID from_point = 0;  //TODO depot
        while(true)
        {
            std::cout << "start routing throw " << cur_core << " core" << std::endl;
            std::cout << "begin_offset is " << begin_offset << std::endl;
            std::cout << "end_offset is " << end_offset << std::endl;
            std::cout << "from_point is " << from_point << std::endl;
            std::cout << "start_point is " << start_point << std::endl;
            Tour tempTour(chains, time_matrix);
            coreCWVisitation(tempTour, cur_core, from_point, start_point);
            std::cout << "CWVisitation is" << std::endl << tempTour << std::endl;
            if(tryInsertRange(tempTour.begin(), tempTour.end(), begin_offset, end_offset))
            {
                std::cout << "start inserting inner core points" << std::endl;
                cores[cur_core].attend();
                for(ChainID point : cores[cur_core].getInners())
                    if(!chains[point]->isAttend())
                        tryInsert(point, begin_offset, end_offset);
                    
                std::cout << "core solution is" << std::endl << tours.back() << std::endl;
                std::cout << "start inserting tails" << std::endl;
                    
                for(ChainID point : cores[cur_core].getTails())
                    if(!chains[point]->isAttend() && 
                        std::all_of(chains[point]->getCores().begin(), 
                                    chains[point]->getCores().end(), 
                                    [&](CoreID c){ return cores[c].isAttend(); }))
                    {
                        unsigned max_lbound=0;
                        unsigned max_rbound=0;
                        for(ChainID gate : chains[point]->getGates())
                        {
                            unsigned lbound=0, rbound=0;
                            tours.back().findBounds(gate, lbound, rbound);
                            
                            if(max_lbound < lbound) max_lbound = lbound;
                            if(max_rbound < rbound) max_rbound = rbound;
                        }
                        tryInsert(point, max_lbound, max_rbound);
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
                        for(ChainID out_gate : chains[cross_tail]->getGates())
                            if(chains[out_gate]->isInner(core))
                                for(ChainID in_gate : chains[cross_tail]->getGates())
                                    if(!cores[chains[in_gate]->getCoreInner()].isAttend())
                                        if(min_dist > chains[out_gate]->timeTo(in_gate))
                                        {
                                            min_dist = chains[out_gate]->timeTo(in_gate);
                                            next_core = chains[in_gate]->getCoreInner();
                                            min_out_gate = out_gate;
                                            min_in_gate = in_gate;
                                            from_core = core;
                                        }
            if(min_dist == UNREACHED_WEIGHT)
                break;
            
            std::cout << "best coice is " << next_core << " conected whis " << from_core << " throw " << min_out_gate << "-" << min_in_gate << std::endl;
            std::cout << "start looking for insert bounds" << std::endl;
            

            unsigned max_lbound=0;
            unsigned max_rbound=0;
            for(ChainID cross_tail : cores[from_core].getTails())
                if(chains[cross_tail]->hasCore(next_core))
                    for(ChainID out_gate : chains[cross_tail]->getGates())
                        if(chains[out_gate]->isInner(from_core))
                        {
                            unsigned lbound=0, rbound=0;
                            tours.back().findBounds(out_gate, lbound, rbound);
                            if(max_lbound < lbound) max_lbound = lbound;
                            if(max_rbound < rbound) max_rbound = rbound;
                        }
            begin_offset = max_lbound;
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
                  const std::vector<FixedPointCoordinate> &coordinates_container)
        : n_points(n_points),
          time_matrix(time_matrix),
          length_matrix(length_matrix),
          full_point_graph(n_points)
    {
        int i=0;
        for(const auto &row : full_forward_graph_container)
            full_point_graph.insert(i++, row.begin(), row.end());
        for(auto &coor : coordinates_container)
            coordinates.emplace_back(coor.lat/1000000.0, coor.lon/1000000.0);
    }
    
    void run()
    {
        findChains();
        buildNearestGraph();
        findCorePoints();
        //disunitCores();
        findCoreForwardTails();
        findCoreBackwardTails();
        linkUnattachedCores();
        
        tours.emplace_back(chains, time_matrix);
        chains[0]->attend();  //TODO depot
        
        CoreID max_core;
        ChainID max_chain;
        int max_dist = 0;
        for(CoreID c=0; c<n_cores; ++c)
        {
            int min_dist = UNREACHED_WEIGHT;
            ChainID min;
            for(ChainID v : cores[c].getInners())
                if(min_dist > chains[0]->timeTo(v))
                {
                    min_dist = chains[0]->timeTo(v);
                    min = v;
                }
            if(max_dist < min_dist)
            {
                max_dist = min_dist;
                max_core = c;
                max_chain = min;
            }
        }
        startRoutingFromCore(max_core, max_chain);
    }
    
    void render(std::vector<char> &output)
    {
        JSON::Object json_root;
        JSON::Array tours_array;
        JSON::Array nearest_array;
        JSON::Array chains_array;
        JSON::Array cores_array;
        for(auto tour : tours)
        {
            JSON::Array tour_points_array;
            tour_points_array.values.push_back(0); //TODO depot
            for(ChainID p : tour)
                std::copy(chains[p]->begin(), chains[p]->end(), std::back_inserter(tour_points_array.values));
            tour_points_array.values.push_back(0); //TODO depot
            tours_array.values.push_back(tour_points_array);
        }
        json_root.values["tours"] = tours_array;
        for(ChainID start=0; start<n_chains; ++start)
        {
            JSON::Array nearest_array_row;
            for(const ChainID point : nearest_graph.forward(start))
                nearest_array_row.values.push_back(point);
            nearest_array.values.push_back(nearest_array_row);
        }
        json_root.values["nearest"] = nearest_array;
        for(ChainID p=0; p<n_chains; ++p)
        {
            JSON::Array chains_array_row;
            std::copy(chains[p]->begin(), chains[p]->end(), std::back_inserter(chains_array_row.values));
            nearest_array.values.push_back(chains_array_row);
        }
        json_root.values["chains"] = chains_array;
        for(CoreID c=0; c<n_cores; ++c)
        {
            JSON::Array cores_array_row;
            std::copy(cores[c].getInners().begin(), cores[c].getInners().end(), std::back_inserter(cores_array_row.values));
            nearest_array.values.push_back(cores_array_row);
        }
        json_root.values["cores"] = cores_array;
        json_root.values["n"] = n_points;
        json_root.values["n_chains"] = n_chains;
        json_root.values["n_cores"] = n_cores;
        JSON::render(output, json_root);
        
    }
};
    
#endif //GRAPG_LOGISTIC_H