//
// Created by maksym on 23.11.16.
//

#include "engine/plugins/math.hpp"
#include "engine/plugins/table.hpp"

#include "engine/api/table_api.hpp"
#include "engine/api/table_parameters.hpp"
#include "engine/api/route_parameters.hpp"
#include "engine/routing_algorithms/many_to_many.hpp"
#include "engine/search_engine_data.hpp"
#include "util/json_container.hpp"
#include "util/string_util.hpp"
#include "engine/guidance/leg_geometry.hpp"
#include "engine/guidance/assemble_geometry.hpp"

#include <cstdlib>

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include <boost/assert.hpp>
#include <engine/api/route_parameters.hpp>

namespace osrm
{
    namespace engine
    {
        namespace plugins
        {

            MathPlugin::MathPlugin(datafacade::BaseDataFacade &facade, const int max_locations_distance_table)
                    : BasePlugin{facade}, distance_table(&facade, heaps), shortest_path(&facade, heaps),
                      max_locations_distance_table(max_locations_distance_table)
            {
            }

            private:
            std::shared_ptr<SearchEngine<DataFacadeT>> search_engine_ptr;


            public:
            explicit MathPlugin(DataFacadeT *facade) : descriptor_string("math"), facade(facade)
                    {
                            search_engine_ptr = std::make_shared<SearchEngine<DataFacadeT>>(facade);
                    }

            Status MathPlugin::HandleRequest(const api::TableParameters &params, util::json::Object &result)
            {
                BOOST_ASSERT(params.IsValid());

                if (!CheckAllCoordinates(params.coordinates))
                {
                    return Error("InvalidOptions", "Coordinates are invalid", result);
                }

                if (params.bearings.size() > 0 && params.coordinates.size() != params.bearings.size())
                {
                    return Error(
                            "InvalidOptions", "Number of bearings does not match number of coordinates", result);
                }

                // Empty sources or destinations means the user wants all of them included, respectively
                // The ManyToMany routing algorithm we dispatch to below already handles this perfectly.
                const auto num_sources =
                        params.sources.empty() ? params.coordinates.size() : params.sources.size();
                const auto num_destinations =
                        params.destinations.empty() ? params.coordinates.size() : params.destinations.size();

                if (max_locations_distance_table > 0 &&
                    ((num_sources * num_destinations) >
                     static_cast<std::size_t>(max_locations_distance_table * max_locations_distance_table)))
                {
                    return Error("TooBig", "Too many table coordinates", result);
                }

                auto snapped_phantoms = SnapPhantomNodes(GetPhantomNodes(params));
                auto result_table = distance_table(snapped_phantoms, params.sources, params.destinations);

                if (result_table.first.empty())
                {
                    return Error("NoTable", "No table found", result);
                }

                api::TableAPI table_api{facade, params};
                table_api.MakeResponse(result_table, snapped_phantoms, result);

                // Check number of parameters
                if (1 > params.coordinates.size())
                {
                    return Error(
                            "InvalidOptions", "Number of parameters is too smal", result);
                }

                RawRouteData raw_route;
                raw_route.check_sum = facade.GetCheckSum();


                if (std::any_of(begin(params.coordinates),
                                end(params.coordinates),
                                [&](Coordinate coordinate)
                                { return !coordinate.IsValid(); }))
                {
                    return Error(
                            "InvalidOptions", "Coordinate is not valid", result);
                }

                for (const Coordinate &coordinate : params.coordinates)
                {
                    raw_route.raw_via_node_coordinates.emplace_back(std::move(coordinate));
                }

                const bool checksum_OK = (params.check_sum == raw_route.check_sum);
                unsigned max_locations = static_cast<unsigned>(raw_route.raw_via_node_coordinates.size());
                PhantomNodeArray phantom_node_vector(max_locations);
                for (unsigned i = 0; i < max_locations; ++i)
                {
                     if (checksum_OK && i < params.hints.size() &&
                       !params.hints[i].empty())
                     {
                         PhantomNode current_phantom_node;
                           DecodeObjectFromBase64(params.hints[i], current_phantom_node);
                        if (current_phantom_node.isValid(facade.GetNumberOfNodes()))
                        {
                            phantom_node_vector[i].emplace_back(std::move(current_phantom_node));
                            continue;
                        }
                    }
                    facade.IncrementalFindPhantomNodeForCoordinate(raw_route.raw_via_node_coordinates[i],
                                                                       phantom_node_vector[i],
                                                                       params.zoom_level,
                                                                       1);


                    BOOST_ASSERT(phantom_node_vector[i].front().isValid(facade->GetNumberOfNodes()));
                }
                search_engine_ptr->math(phantom_node_vector,
                                           params.coordinates,
                                           params.transport_restrictions.front(),
                                           reply.content);

                return Status::Ok;
            }
        }
    }
}
