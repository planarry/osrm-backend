#ifndef GRAPG_H
#define GRAPG_H

#include "CommonHeader.h"

template<typename T>
class Graph {
    unsigned n;
    std::vector<std::set<T>> forward_data, reverse_data;
    std::vector<std::list<T>> forward_order_data;
    
public:
    Graph() : n(0) { }
    Graph(unsigned n) : n(n), forward_data(n), reverse_data(n), forward_order_data(n) { }
    
    void assign(unsigned n)
    {
        this->n = n;
        forward_data.assign(n, std::set<T>());
        reverse_data.assign(n, std::set<T>());
        forward_order_data.assign(n, std::list<T>());
    }
    
    void insert(T from, T to)
    {
        BOOST_ASSERT(from < n);
        BOOST_ASSERT(to < n);
        if(!has(from, to)) //TODO depot
        {
            forward_data[from].insert(to);
            reverse_data[to].insert(from);
            if(to != 0)
                forward_order_data[from].push_back(to);
        }
    }
    
    template<typename Iterator>
    void insert(T from, Iterator begin, Iterator end)
    {
        BOOST_ASSERT(from < n);
        std::for_each(begin, end, [&] (T to) {
            insert(from, to);
        });
    }
    
    void erase(T from, T to)
    {
        BOOST_ASSERT(from < n);
        BOOST_ASSERT(to < n);
        if(!has(from, to))
        {
            forward_data[from].erase(to);
            reverse_data[to].erase(from);
            forward_order_data[from].erase(std::find(forward_order_data[from].begin(), 
                                                     forward_order_data[from].end(),
                                                     to));
        }
    }
        
    bool has(T from, T to) const
    { 
        BOOST_ASSERT(from < n);
        BOOST_ASSERT(to < n);
        return forward_data[from].find(to) != forward_data[from].end();
    }
    
    const std::list<T>& forward(T from) const
    { 
        BOOST_ASSERT(from < n);
        return forward_order_data[from]; 
    }
    
    const std::set<T>& forward_set(T from) const
    {
        BOOST_ASSERT(from < n);
        return forward_data[from]; 
    }
    
    const std::set<T>& reverse(T from) const
    {
        BOOST_ASSERT(from < n);
        return reverse_data[from]; 
    }
};

#endif