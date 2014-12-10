#ifndef CORE_H
#define CORE_H

#include "CommonHeader.h"

class Core {
    CoreID id;
    bool is_attend;
    std::set<ChainID> inners_data, tails_data;
    std::vector<std::shared_ptr<ChainBase>> &chains;
    std::vector<Core> &cores;
public:
    Core(std::vector<std::shared_ptr<ChainBase>> &chains, std::vector<Core> &cores)
        : is_attend(false), chains(chains), cores(cores)
    {}
    
    void insert(ChainID inner)
    { inners_data.insert(inner); }
    
    void insertTail(ChainID tail)
    { tails_data.insert(tail); }
    
    unsigned size() const
    { return inners_data.size(); }
    
    const std::set<ChainID>& getInners() const
    { return inners_data; }
    
    const std::set<ChainID>& getTails() const
    { return tails_data; }
    
    void attend()
    { is_attend = true; }
    
    bool isAttend()
    { return is_attend; }
};

#endif