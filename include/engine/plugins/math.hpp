//
// Created by maksym on 23.11.16.
//

#ifndef OSRM_MATH_HPP
#define OSRM_MATH_HPP


#include <engine/api/route_parameters.hpp>
#include "engine/plugins/plugin_base.hpp"

#include "engine/api/table_parameters.hpp"
#include "engine/routing_algorithms/many_to_many.hpp"
#include "engine/routing_algorithms/shortest_path.hpp"
#include "engine/search_engine_data.hpp"
#include "util/json_container.hpp"

#include "boost/optional.hpp"

namespace osrm
{
namespace engine
{
namespace plugins
{

class MathPlugin final : public BasePlugin
{
  public:
    explicit MathPlugin(datafacade::BaseDataFacade &facade,
                         const int max_locations_distance_table);

    Status HandleRequest(const api::TableParameters &params, util::json::Object &result);

  private:
    SearchEngineData heaps;
    routing_algorithms::ManyToManyRouting<datafacade::BaseDataFacade> distance_table;
    routing_algorithms::ShortestPathRouting<datafacade::BaseDataFacade> shortest_path;
    int max_locations_distance_table;
    /*std::shared_ptr<SearchEngine<DataFacadeT>> search_engine_ptr;
    std::string descriptor_string;
    DataFacadeT *facade;*/
};
}
}
}


#endif //OSRM_MATH_HPP
