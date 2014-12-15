#ifndef CHAIN_H
#define CHAIN_H

#include "CommonHeader.h"
#include <boost/iterator/iterator_facade.hpp>
#include <boost/variant/variant.hpp>


class chain_iterator
  : public boost::iterator_facade< chain_iterator , const PointID , boost::random_access_traversal_tag>
{
public:
    typedef boost::variant<const PointID*, std::vector<PointID>::const_iterator, std::vector<PointID>::const_reverse_iterator> base_iterator;
    
    chain_iterator() {}

    chain_iterator(const PointID* base)
      : m_base(base) {}

    chain_iterator(std::vector<PointID>::const_iterator base)
      : m_base(base) {}

    chain_iterator(std::vector<PointID>::const_reverse_iterator base)
      : m_base(base) {}

private:
    friend class boost::iterator_core_access;

    struct do_increment: boost::static_visitor<void>
    {
        template<typename T>
        void operator()( T& t ) const
        { ++t; }
    };

    struct do_decrement: boost::static_visitor<void>
    {
        template<typename T>
        void operator()( T& t ) const
        { --t; }
    };

    struct do_advance: boost::static_visitor<void>
    {
        int _n;
        do_advance(int n) : _n(n) {}
        template<typename T>
        void operator()( T& t) const
        { t+=_n; }
    };

    struct do_distance_to : boost::static_visitor<unsigned>
    {
        template <typename A, typename B>
        unsigned operator() ( const A & a, const B & b) const
        { return 0; }

        template <typename A>
        unsigned operator() ( const A & a, const A & b) const
        { return b-a; }
    };

    struct do_compare : boost::static_visitor<bool>
    {
        template <typename A, typename B>
        bool operator() ( const A & a, const B & b) const
        { return true; }

        template <typename A>
        bool operator() ( const A & a, const A & b) const
        { return a == b; }
    };
    
    struct do_dereference: boost::static_visitor<const PointID&>
    {
        template<typename T>
        const PointID& operator()( const T& t ) const
        { return *t; }
    };
    
    void increment() 
    { boost::apply_visitor(do_increment(), m_base); }
    
    void decrement() 
    { boost::apply_visitor(do_decrement(), m_base); }
    
    void advance(int n) 
    { boost::apply_visitor(do_advance(n), m_base); }
    
    unsigned distance_to(const chain_iterator & other) const
    { return boost::apply_visitor(do_distance_to(), m_base, other.m_base); }
    
    bool equal(const chain_iterator & other) const
    { return boost::apply_visitor(do_compare(), m_base, other.m_base); }

    const PointID& dereference() const
    { return boost::apply_visitor(do_dereference(), m_base); }

    base_iterator m_base;
};

class ReverseChain;
class ChainBase 
{
protected:
    friend ReverseChain;
    ChainID id;
    CoreID coreInner;
    bool is_attend;
    int inner_time;
    std::set<CoreID> cores_data;
    std::set<ChainID> gates_data;
    std::vector<std::shared_ptr<ChainBase>> &chains;
    std::vector<Core> &cores;
    const ublas::matrix<int> &time_matrix, &length_matrix;
    
public:    
    virtual chain_iterator begin() const = 0;
    
    virtual chain_iterator end() const = 0;
    
    ChainBase(ChainID id,
              std::vector<std::shared_ptr<ChainBase>> &chains,
              std::vector<Core> &cores,
              const ublas::matrix<int> &time_matrix,
              const ublas::matrix<int> &length_matrix)
        : id(id),
          coreInner(SPECIAL_ID),
          is_attend(false), 
          inner_time(UNREACHED_WEIGHT),
          chains(chains),
          cores(cores),
          time_matrix(time_matrix),
          length_matrix(length_matrix)
    { }
    
    void setInner(CoreID core)
    {
        coreInner = core;
        cores_data.insert(core);
        gates_data.insert(id);
    }
    CoreID getCoreInner()
    { return coreInner; }
    
    bool hasCore(CoreID core)
    { return cores_data.find(core) != cores_data.end(); }
    
    bool isUnattached()
    { return cores_data.empty(); }
    
    const std::set<CoreID>& getCores() const
    { return cores_data; }
    
    const std::set<CoreID>& getGates() const
    { return gates_data; }
    
    void bindTail(ChainID to)
    {
        for(const CoreID &c : chains[to]->getCores())
        {
            cores[c].insertTail(id);
            cores_data.insert(c);
        }
        gates_data.insert(chains[to]->gates_data.begin(), chains[to]->gates_data.end());
    }
    
    void bindTail(CoreID core, ChainID gate)
    {
        cores[core].insertTail(id);
        cores_data.insert(core);
        gates_data.insert(gate);
    }
    
    int timeTo(ChainID to)
    {
        int time = time_matrix(this->getOutPointID(), chains[to]->getInPointID());
        if(time < time_matrix(this->getOutPointID(), chains[chains[to]->reverse()]->getInPointID()))
            --time;
        return time;
    }
    
    bool isInner() const
    { return coreInner != SPECIAL_ID; }
    
    bool isInner(CoreID core) const
    { return coreInner == core; }
    
    bool sameCoreInner(ChainID with) const
    { return isInner() && coreInner==chains[with]->coreInner; }
    
    int getInnerTime()
    { 
        if(inner_time == UNREACHED_WEIGHT)
        {
            inner_time = 0;
            for(auto i = ++this->begin(); i<this->end(); ++i)
            {
                inner_time += time_matrix(*(i - 1), *i);
                for(unsigned j=1; j<=MOVING_MEAN_RADIUS; ++j)
                {
                    if(i - j > this->begin())
                        inner_time += time_matrix(*(i - j - 1), *(i)) / (j * MOVING_MEAN_FACTOR);
                    if(i + j < this->end())
                        inner_time += time_matrix(*(i - 1), *(i + j)) / (j * MOVING_MEAN_FACTOR);
                }
            }
        }
    
        return inner_time; 
    }
    
    bool isAttend() const
    { return is_attend; }
    
    virtual void attend()
    { is_attend = true; }
    
    virtual ChainID reverse() const
    { return id; }
    
    virtual PointID getInPointID() const = 0;
    
    virtual PointID getOutPointID() const = 0;
};

class Point : public ChainBase {
    PointID point;
public:
    Point(ChainID id,
          PointID point,
          std::vector<std::shared_ptr<ChainBase>> &chains,
          std::vector<Core> &cores,
          const ublas::matrix<int> &time_matrix,
          const ublas::matrix<int> &length_matrix)
        : ChainBase(id, chains, cores, time_matrix, length_matrix), point(point) 
    { }
    
    virtual PointID getInPointID() const
    { return point; }
    
    virtual PointID getOutPointID() const
    { return point; }
    
    virtual chain_iterator begin() const
    { return &point; }
    
    virtual chain_iterator end() const
    { return &point + 1; }
};

class OnedirChain : public ChainBase {
    std::vector<PointID> points;
public:
    OnedirChain(ChainID id, 
                const std::list<PointID> &_points, 
                std::vector<std::shared_ptr<ChainBase>> &chains, 
                std::vector<Core> &cores,
                const ublas::matrix<int> &time_matrix,
                const ublas::matrix<int> &length_matrix) 
        : ChainBase(id, chains, cores, time_matrix, length_matrix), points(_points.begin(), _points.end()) 
    { }
    
    virtual PointID getInPointID() const
    { return points.front(); }
    
    virtual PointID getOutPointID() const
    { return points.back(); }
    
    virtual chain_iterator begin() const
    { return points.cbegin(); }
    
    virtual chain_iterator end() const
    { return points.cend(); }
};

class ForwardChain : public ChainBase {
    friend ReverseChain;
    ChainID reverse_id;
    std::vector<PointID> points;
public:
    ForwardChain(ChainID id, 
                 const std::list<PointID> &_points, 
                 std::vector<std::shared_ptr<ChainBase>> &chains, 
                 std::vector<Core> &cores,
                 const ublas::matrix<int> &time_matrix,
                 const ublas::matrix<int> &length_matrix) 
        : ChainBase(id, chains, cores, time_matrix, length_matrix), reverse_id(SPECIAL_ID), points(_points.begin(), _points.end())
    { }
    
    virtual PointID getInPointID() const
    { return points.front(); }
    
    virtual PointID getOutPointID() const
    { return points.back(); }
    
    virtual void attend()
    { 
        BOOST_ASSERT(reverse_id != SPECIAL_ID);
        chains[reverse_id]->attend(); 
    }
    
    virtual ChainID reverse() const
    { return reverse_id; }
    
    virtual chain_iterator begin() const
    { return points.cbegin(); }
    
    virtual chain_iterator end() const
    { return points.cend(); }
};

class ReverseChain : public ChainBase {
    ForwardChain &base;
public:
    ReverseChain(ChainID id, ForwardChain &forwardChain)
        : ChainBase(id, forwardChain.chains, forwardChain.cores, forwardChain.time_matrix, forwardChain.length_matrix), 
          base(forwardChain)
    { forwardChain.reverse_id = id; }
    
    virtual PointID getInPointID() const
    { return base.getOutPointID(); }
    
    virtual PointID getOutPointID() const
    { return base.getInPointID(); }
    
    virtual void attend()
    { 
        is_attend = true; 
        base.is_attend = true;
    }
    
    virtual ChainID reverse() const
    { return base.id; }
    
    virtual chain_iterator begin() const
    { return base.points.crbegin(); }
    
    virtual chain_iterator end() const
    { return base.points.crend(); }
};

#endif