
#ifndef TRANSPORT_RESTRICTION_H
#define TRANSPORT_RESTRICTION_H

#include "../../Util/SimpleLogger.h"
#include "../../typedefs.h"

#include <vector>
#include <boost/spirit/include/qi.hpp>

struct TransportRestriction
{
  short load;
  short height;
  
  TransportRestriction(short load, short height):load(load),height(height){
    //SimpleLogger().Write()<<"load="<<load;
    //SimpleLogger().Write()<<"height="<<height;
  }
  TransportRestriction():load(0),height(0){}
  
  template <typename EdgeDataT>
  bool IsEdgeRestricted(EdgeDataT data) const
  {
    //SimpleLogger().Write()<<"load="<<load;
    //SimpleLogger().Write()<<"height="<<height;
    //SimpleLogger().Write()<<"maxload="<<data.maxload;
    //SimpleLogger().Write()<<"maxheight="<<data.maxheight;
    return (data.maxload<load) || (data.maxheight<height);
  }
  
  template <typename Iterator>
  static boost::spirit::qi::rule<Iterator> getRule()
  {
    return boost::spirit::qi::double_ >> boost::spirit::qi::lit(',') >> boost::spirit::qi::double_;
  }
};

#endif // TRANSPORT_RESTRICTION_H
