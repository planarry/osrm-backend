#ifndef TOUR_H
#define TOUR_H

#include "CommonHeader.h"

class Tour {
    std::vector<ChainID> data;
    std::map<ChainID, unsigned> indexes;
    std::vector<std::shared_ptr<ChainBase>> &chains;
    const ublas::matrix<int> &time_matrix;
    int _cost;
public:
    Tour(std::vector<std::shared_ptr<ChainBase>> &chains, const ublas::matrix<int> &time_matrix) 
        : chains(chains), time_matrix(time_matrix), _cost(UNREACHED_WEIGHT)
    {}
    Tour(Tour &&r) 
        : chains(r.chains), time_matrix(r.time_matrix), _cost(UNREACHED_WEIGHT)
    { *this = std::move(r); }
    
    Tour(Tour &r) 
        : data(r.data), indexes(r.indexes), chains(r.chains), time_matrix(r.time_matrix), _cost(r._cost)
    {}
    
    Tour(const Tour &templateTour, int add_count)
        : chains(templateTour.chains), time_matrix(templateTour.time_matrix), _cost(UNREACHED_WEIGHT)
    { data.reserve(templateTour.data.size() + add_count); }
    
    typedef typename std::vector<ChainID>::iterator iterator;
    typedef typename std::vector<ChainID>::const_iterator const_iterator;
    
    const_iterator begin(unsigned lbound = 0) const
    { return data.cbegin() + lbound; }
    const_iterator end(unsigned rbound = 0) const
    { return data.cend() - rbound; }
    iterator begin(unsigned lbound = 0)
    { return data.begin() + lbound; }
    iterator end(unsigned rbound = 0)
    { return data.end() - rbound; }
    unsigned size()
    { return data.size(); }
    
    bool operator>(Tour &r)
    { return cost() > r.cost(); }
    bool operator>=(Tour &r)
    { return cost() >= r.cost(); }
    bool operator<(Tour &r)
    { return cost() < r.cost(); }
    bool operator==(Tour &r)
    { return cost() == r.cost(); }
    
    Tour &operator=(Tour &&r)
    {
        if (this != &r)
        {
            data = std::move(r.data);
            indexes = std::move(r.indexes);
            _cost = r._cost;
        }
        return *this;
    }
    
    int cost(){
        if(_cost == UNREACHED_WEIGHT && !data.empty())
        {
            _cost = chains[0]->timeTo(data.front()) + chains[data.back()]->timeTo(0); //TODO depot
            _cost += chains[data.front()]->getInnerTime();
            for(unsigned i=1; i<data.size(); ++i)
            {
                _cost += chains[data[i]]->getInnerTime();
                _cost += chains[data[i - 1]]->timeTo(data[i]);
                for(unsigned j=1; j<=MOVING_MEAN_RADIUS; ++j)
                {
                    if(int(i - j) > 0)
                        _cost += chains[data[i - j - 1]]->timeTo(data[i]) / (j * MOVING_MEAN_FACTOR);
                    if(i + j < data.size())
                        _cost += chains[data[i - 1]]->timeTo(data[i + j]) / (j * MOVING_MEAN_FACTOR);
                }
            }
        }
        return _cost;
    }
    
    template<typename ForwardIterator>
    void copy(ForwardIterator begin, ForwardIterator end)
    {
        std::for_each(begin, end, [&](ChainID chain){
            this->add(chain);
        });
    }
    template<typename BidirIterator>
    void reverse_copy(BidirIterator begin, BidirIterator end)
    {
        std::for_each(std::reverse_iterator<BidirIterator>(end), std::reverse_iterator<BidirIterator>(begin), [&](ChainID chain){
            this->add(chains[chain]->reverse());
        });
    }
    void findBounds(ChainID gate, unsigned int &lbound, unsigned int &rbound)
    {
        auto is_core_point=[&](ChainID i) { return chains[i]->isInner(); };
        auto gate_iter = begin(indexes[gate]);
        rbound = end() - std::find_if(gate_iter + 1, end(), is_core_point);
        //if(rbound != tours.back().end() - end_offset) ++rbound;
        lbound = std::find_if(std::reverse_iterator<iterator>(gate_iter-1), 
                                std::reverse_iterator<iterator>(begin()), 
                                is_core_point).base() - begin();
    }
    void add(ChainID chain)
    {
        indexes[chain] = data.size();
        data.push_back(chain);
    }
};

std::ostream& operator<< (std::ostream& out, const Tour &r) {
    std::copy(r.begin(), r.end(), std::ostream_iterator<ChainID>(out, " "));
    return out;
}

#endif