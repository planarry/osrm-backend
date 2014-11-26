#ifndef MATH_PLUGIN_H
#define MATH_PLUGIN_H

#include "BasePlugin.h"

#include "../Algorithms/ObjectToBase64.h"
#include "../DataStructures/JSONContainer.h"
#include "../DataStructures/QueryEdge.h"
#include "../DataStructures/SearchEngine.h"
#include "../Descriptors/BaseDescriptor.h"
#include "../Util/SimpleLogger.h"
#include "../Util/StringUtil.h"
#include "../Util/TimingUtil.h"

#include <cstdlib>

#include <algorithm>
#include <memory>
#include <unordered_map>
#include <string>
#include <vector>

template <class DataFacadeT> class MathPlugin : public BasePlugin
{
  private:
    std::shared_ptr<SearchEngine<DataFacadeT>> search_engine_ptr;

  public:
    explicit MathPlugin(DataFacadeT *facade) : descriptor_string("math"), facade(facade)
    {
        search_engine_ptr = std::make_shared<SearchEngine<DataFacadeT>>(facade);
    }

    virtual ~MathPlugin() {}

    const std::string GetDescriptor() const { return descriptor_string; }

    void HandleRequest(const RouteParameters &route_parameters, http::Reply &reply)
    {
        TIMER_START(request);
        // check number of parameters
        if (1 > route_parameters.coordinates.size())
        {
            reply = http::Reply::StockReply(http::Reply::badRequest);
            return;
        }

        RawRouteData raw_route;
        raw_route.check_sum = facade->GetCheckSum();

        if (std::any_of(begin(route_parameters.coordinates),
                        end(route_parameters.coordinates),
                        [&](FixedPointCoordinate coordinate)
                        { return !coordinate.isValid(); }))
        {
            reply = http::Reply::StockReply(http::Reply::badRequest);
            return;
        }

        for (const FixedPointCoordinate &coordinate : route_parameters.coordinates)
        {
            raw_route.raw_via_node_coordinates.emplace_back(std::move(coordinate));
        }

        const bool checksum_OK = (route_parameters.check_sum == raw_route.check_sum);
        unsigned max_locations = static_cast<unsigned>(raw_route.raw_via_node_coordinates.size());
        PhantomNodeArray phantom_node_vector(max_locations);
        for (unsigned i = 0; i < max_locations; ++i)
        {
            if (checksum_OK && i < route_parameters.hints.size() &&
                !route_parameters.hints[i].empty())
            {
                PhantomNode current_phantom_node;
                DecodeObjectFromBase64(route_parameters.hints[i], current_phantom_node);
                if (current_phantom_node.isValid(facade->GetNumberOfNodes()))
                {
                    phantom_node_vector[i].emplace_back(std::move(current_phantom_node));
                    continue;
                }
            }
            facade->IncrementalFindPhantomNodeForCoordinate(raw_route.raw_via_node_coordinates[i],
                                                            phantom_node_vector[i],
                                                            route_parameters.zoom_level,
                                                            1);

            BOOST_ASSERT(phantom_node_vector[i].front().isValid(facade->GetNumberOfNodes()));
        }
        search_engine_ptr->math(phantom_node_vector, 
                                route_parameters.coordinates,
                                route_parameters.transport_restrictions.front(), 
                                reply.content);
        //for(auto p : *points)
        //    SimpleLogger().Write()<<p;
        TIMER_STOP(request);
        SimpleLogger().Write() << "Request processing time is " << TIMER_SEC(request) << "s";
    }

  private:
    std::string descriptor_string;
    DataFacadeT *facade;
};

#endif // MATH_PLUGIN_H
