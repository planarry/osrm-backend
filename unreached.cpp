//
// Created by dan1992 on 20.11.15.
//

#include "DataStructures/ImportEdge.h"
#include "DataStructures/ImportNode.h"
#include "Util/GraphLoader.h"
#include "Util/FingerPrint.h"
#include "typedefs.h"

#include <boost/filesystem/fstream.hpp>
#include <boost/program_options.hpp>
#include <string>
#include <vector>
#include <stack>
#include <map>
#include <set>


std::string restrictions_path = "map.osrm.restrictions";
std::string path = "map.osrm";

std::string out_restrictions_path = "new_map.osrm.restrictions";
std::string out_path = "new_map.osrm";

std::vector<NodeInfo> internal_to_external_node_map;
std::vector<TurnRestriction> restriction_list;
std::vector<NodeID> barrier_node_list;
std::vector<NodeID> traffic_light_list;
std::vector<ImportEdge> edge_list;

std::vector<bool> forward_visited, backward_visited, used_nodes;
std::vector<std::set<unsigned int>> nodes;
std::map<NodeID, std::set<unsigned int>> restrictions;

struct Config {
    std::string inputFile;
    std::string outputFile;
    NodeID startNode;
};

bool ParseArguments(int argc, char *argv[], Config& config) {
    boost::program_options::options_description config_options("Configuration");
    config_options.add_options()(
            "inputFile,i",
            boost::program_options::value<std::string>(&config.inputFile),
            "input *.osrm filename"
    )(
            "outputFile,o",
            boost::program_options::value<std::string>(&config.outputFile),
            "output *.osrm filename"
    )(
            "startNode,n",
            boost::program_options::value<unsigned int>(&config.startNode)->default_value(1465728853),
            "ID of starting node"
    );                                                                                  //3618470742 - test 1465728853

    boost::program_options::variables_map option_variables;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv)
                                          .options(config_options)
                                          .run(),
                                  option_variables);

    boost::program_options::notify(option_variables);

    if (!option_variables.count("inputFile")) {
        std::cout << "missing input file" << std::endl;
        return false;
    }

    if (!option_variables.count("outputFile")) {
        config.outputFile = config.inputFile;
    }

    return true;
}

void CheckRestrictionsFile(FingerPrint &fingerprint_orig) {
    boost::filesystem::ifstream restriction_stream(restrictions_path, std::ios::binary);
    FingerPrint fingerprint_loaded;
    unsigned number_of_usable_restrictions = 0;
    restriction_stream.read((char *) &fingerprint_loaded, sizeof(FingerPrint));
    if (!fingerprint_loaded.TestPrepare(fingerprint_orig)) {
        SimpleLogger().Write(logWARNING) << ".restrictions was prepared with different build.\n"
                "Reprocess to get rid of this warning.";
    }

    restriction_stream.read((char *) &number_of_usable_restrictions, sizeof(unsigned));
    restriction_list.resize(number_of_usable_restrictions);
    if (number_of_usable_restrictions > 0) {
        restriction_stream.read((char *) &(restriction_list[0]),
                                number_of_usable_restrictions * sizeof(TurnRestriction));
    }
    restriction_stream.close();
}

bool checkTurnRestriction(const NodeID from, const NodeID via, const NodeID to) {
    if (restrictions.find(via) == restrictions.end())
        return true;
    bool allow = true;
    TurnRestriction restriction;
    for (unsigned int i: restrictions[via]) {
        restriction = restriction_list[i];
        if ((restriction.fromNode == from) && (restriction.toNode == to))
            return restriction.flags.isOnly;
        if (restriction.flags.isOnly && (restriction.fromNode == from))
            allow = false;
    }
    return allow;
}

NodeID findNextEdge(NodeID from, NodeID via, const std::vector<bool> &visited, bool backward = false) {
    bool next_source, direction;
//    if((internal_to_external_node_map[via].node_id == 26261475)
//            || (internal_to_external_node_map[via].node_id == 359985524))
//    {
//        std::cout<<"AAAAAAAAAA!";
//    }
    NodeID to, from_t, to_t;
    if (!std::binary_search(barrier_node_list.begin(), barrier_node_list.end(), via)) {     //checking for barrier
        for (unsigned int i: nodes[via]) {                                                  //iterating node edges
            if (visited[i]) continue;                                                       //if visited
            if (edge_list[i].access_restricted) continue;                                   //if restricted
            next_source = (via == edge_list[i].source);                       //checking mutual point
            if (backward) {
                direction = (!next_source) ? edge_list[i].forward : edge_list[i].backward;       //checking direction
            } else {
                direction = (next_source) ? edge_list[i].forward : edge_list[i].backward;
            }
            if (!direction) continue;
            to = next_source ? edge_list[i].target : edge_list[i].source;                   //to node
            if (backward) {
                to_t = from;
                from_t = to;
            } else {
                to_t = to;
                from_t = from;
            }
            if (!checkTurnRestriction(from_t, via, to_t))
                continue;
            return i;
        }
    }
    return UINT_MAX;
}

void depthFirstSearch(const NodeID first, std::vector<bool> &visited, bool backward = false) {
    NodeID current;
    NodeID next, from, via;
    std::stack<NodeID> stack;
    stack.push(first);
    while (!stack.empty()) {                            //forward access
        current = stack.top();

        visited[current] = true;                        //access flag

        if (backward) {
            from = edge_list[current].target;
            via = edge_list[current].source;
        } else {
            from = edge_list[current].source;
            via = edge_list[current].target;
        }

        if (edge_list[current].forward) {               //checking current edge direction
            next = findNextEdge(from, via, visited, backward);
            if (next != UINT_MAX) {
                stack.push(next);
                continue;
            }
        }

        if (edge_list[current].backward) {              //checking current edge direction
            next = findNextEdge(via, from, visited, backward);
            if (next != UINT_MAX) {
                stack.push(next);
                continue;
            }
        }

        stack.pop();
    }
}

int main(int argc, char *argv[]) {
    Config config;
    if (!ParseArguments(argc, argv, config))
        return 0;

    path = config.inputFile + ".osrm";
    restrictions_path = path + ".restrictions";
    out_path = config.outputFile + ".osrm";
    out_restrictions_path = out_path + ".restrictions";

    FingerPrint fingerPrint;
    CheckRestrictionsFile(fingerPrint);

    boost::filesystem::ifstream input_stream(path, std::ios::binary);
    NodeID n = readBinaryOSRMGraphFromStream(input_stream,
                                             edge_list,
                                             barrier_node_list,
                                             traffic_light_list,
                                             &internal_to_external_node_map,
                                             restriction_list);
    nodes.resize(n);
    for (unsigned int i = 0; i < edge_list.size(); ++i) {
        nodes[edge_list[i].source].insert(i);
        nodes[edge_list[i].target].insert(i);
    }
    for (unsigned int i = 0; i < restriction_list.size(); ++i) {
        restrictions[restriction_list[i].viaNode].insert(i);
    }

//    start node ID 3618470742
//    start edge 14320

//    std::cout << n << ' ' << edge_list.size() << ' ' << barrier_node_list.size() << ' ' << traffic_light_list.size() <<
//    ' ' << internal_to_external_node_map.size() <<
//    ' ' << restriction_list.size() << std::endl;

//    for (unsigned int i = 0; i < internal_to_external_node_map.size(); ++i) {
//        std::cout << internal_to_external_node_map[i].node_id << ' ' <<
//        internal_to_external_node_map[i].lat / 1000000.0 << ' ' <<
//        internal_to_external_node_map[i].lon / 1000000.0 << std::endl;
//        for (NodeID j:nodes[i]) {
//            std::cout << edge_list[j].name_id << ' ' <<
//            edge_list[j].source << ' ' <<
//            edge_list[j].target << std::endl;
//        }
//    }

//    for (NodeID i: nodes[n - 2]) {
//        std::cout << i << std::endl;
//    }

//    std::cout << nodes.size() << std::endl <<
//            sizeof(FingerPrint) + 2 * sizeof(unsigned) + internal_to_external_node_map.size() * sizeof(NodeInfo) +
//            edge_list.size() * 33 << std::endl <<
//            sizeof(ImportEdge) << '=' << 3 * sizeof(unsigned) + 2 * sizeof(int) + 4 * sizeof(short) + 5 * sizeof(bool);

    forward_visited.resize(edge_list.size(), false);
    backward_visited.resize(edge_list.size(), false);

    NodeID start_edge = UINT_MAX;
    for (unsigned int i = 0; i < internal_to_external_node_map.size(); i++) {
        if (internal_to_external_node_map[i].node_id == config.startNode) {
            for (NodeID edge: nodes[i]) {
                start_edge = edge;
                break;
            }
            break;
        }
    }
    if (start_edge == UINT_MAX) {
        std::cout << "error: edge not found, try another starting node" << std::endl;
        return 0;
    }
    std::cout << "start edge ID = " << start_edge << std::endl;

    depthFirstSearch(start_edge, forward_visited);           //forward access DFS

    depthFirstSearch(start_edge, backward_visited, true);    //backward access DFS

    int m = 0;
    for (unsigned int i = 0; i < forward_visited.size(); i++)
        if (forward_visited[i]) m++;
    std::cout << forward_visited.size() << ' ' << m << ' ';
    m = 0;
    for (unsigned int i = 0; i < backward_visited.size(); i++)
        if (backward_visited[i]) m++;
    std::cout << m << ' ';
    unsigned int new_nodes = 0;
    unsigned int new_edges = 0;
    used_nodes.resize(n, false);
    for (unsigned int i = 0; i < backward_visited.size(); i++)
        if (forward_visited[i] &&
            backward_visited[i]) {
            new_edges++;
            used_nodes[edge_list[i].source] = true;
            used_nodes[edge_list[i].target] = true;
        }
    std::cout << new_edges << std::endl;

    for (unsigned int i = 0; i < n; ++i) {
        if (used_nodes[i]) new_nodes++;
    }
    std::cout << n << ' ' << new_nodes << std::endl;

    std::ofstream out_restrictions(out_restrictions_path, std::ios::binary);
    out_restrictions.write((char *) &fingerPrint, sizeof(FingerPrint));
    out_restrictions.write((char *) &n, sizeof(unsigned));

    unsigned int used_restrictions = 0;

    for (TurnRestriction restriction : restriction_list) {
        if (used_nodes[restriction.viaNode] &&
            used_nodes[restriction.fromNode] &&
            used_nodes[restriction.toNode]) {
            restriction.viaNode = internal_to_external_node_map[restriction.viaNode].node_id;
            restriction.fromNode = internal_to_external_node_map[restriction.fromNode].node_id;
            restriction.toNode = internal_to_external_node_map[restriction.toNode].node_id;
            out_restrictions.write((char *) &(restriction), sizeof(TurnRestriction));
            used_restrictions++;
        }
    }

    out_restrictions.seekp(std::ios::beg + sizeof(FingerPrint));
    out_restrictions.write((char *) &used_restrictions, sizeof(unsigned));

    out_restrictions.close();

    std::cout << used_restrictions << ' ' << sizeof(FingerPrint) + sizeof(unsigned) + used_restrictions * sizeof(TurnRestriction) << std::endl;

    std::ofstream out(out_path, std::ios::binary);
    out.write((char *) &fingerPrint, sizeof(FingerPrint));
    out.write((char *) &n, sizeof(unsigned));
    ExternalMemoryNode node;
    new_nodes = 0;
    for (unsigned int i = 0; i < internal_to_external_node_map.size(); i++) {
        if (used_nodes[i]) {
            node.node_id = internal_to_external_node_map[i].node_id;
            node.lat = internal_to_external_node_map[i].lat;
            node.lon = internal_to_external_node_map[i].lon;
            node.bollard = std::binary_search(barrier_node_list.begin(), barrier_node_list.end(), i);
            node.trafficLight = std::binary_search(traffic_light_list.begin(), traffic_light_list.end(), i);
            new_nodes ++;
            out.write((char *) &node, sizeof(ExternalMemoryNode));
        }
    }
    std::ios::pos_type previous = out.tellp();
    out.seekp(std::ios::beg + sizeof(FingerPrint));
    out.write((char *) &new_nodes, sizeof(unsigned));
    out.seekp(previous);

    out.write((char *) &new_edges, sizeof(unsigned));
    short direction;
    bool temp;
    for (unsigned int i = 0; i < edge_list.size(); i++) {
        if (forward_visited[i] &&
            backward_visited[i]) {
            out.write((char *) &internal_to_external_node_map[edge_list[i].source].node_id, sizeof(unsigned));
            out.write((char *) &internal_to_external_node_map[edge_list[i].target].node_id, sizeof(unsigned));
            out.write((char *) &edge_list[i].length, sizeof(int));
            direction = (short)(edge_list[i].forward ? edge_list[i].forward != edge_list[i].backward : 2);
            out.write((char *) &direction, sizeof(short));
            out.write((char *)&edge_list[i].weight, sizeof(int));
            out.write((char *)&edge_list[i].type, sizeof(short));
            out.write((char *)&edge_list[i].name_id, sizeof(unsigned));
            temp = edge_list[i].roundabout;
            out.write((char *)&temp, sizeof(bool));
            temp = edge_list[i].in_tiny_cc;
            out.write((char *)&temp, sizeof(bool));
            temp = edge_list[i].access_restricted;
            out.write((char *)&temp, sizeof(bool));
            temp = edge_list[i].contra_flow;
            out.write((char *)&temp, sizeof(bool));
            temp = edge_list[i].is_split;
            out.write((char *)&temp, sizeof(bool));
            out.write((char *)&edge_list[i].maxload, sizeof(short));
            out.write((char *)&edge_list[i].maxheight, sizeof(short));
        }
    }

    out.close();

    std::cout << sizeof(FingerPrint) + (2 * sizeof(unsigned)) + (new_nodes * sizeof(ExternalMemoryNode)) +
                         (new_edges * ((5 * sizeof(bool)) + (4 * sizeof(short)) + (3 * sizeof(unsigned)) + (2 * sizeof(int)))) << std::endl <<
            sizeof(FingerPrint) + (2 * sizeof(unsigned)) + (nodes.size() * sizeof(ExternalMemoryNode)) +
            (edge_list.size() * ((5 * sizeof(bool)) + (4 * sizeof(short)) + (3 * sizeof(unsigned)) + (2 * sizeof(int)))) << std::endl <<
            (5 * sizeof(bool)) + (4 * sizeof(short)) + (3 * sizeof(unsigned)) + (2 * sizeof(int)) << std::endl;

//    TurnRestriction restr;
//
//    for (unsigned int i = 0; i < internal_to_external_node_map.size(); i++) {
//        switch (internal_to_external_node_map[i].node_id) {
//            case 769:
//            case 773:
//            case 771:
//            case 765:
//            case 766:
//            case 767:
//            case 777:
//            case 1657834148:
//                if (!used_nodes[i])
//                    std::cout << "+ unaccessable " << internal_to_external_node_map[i].node_id << std::endl;
//                else
//                    std::cout << "- accessable !!! " << internal_to_external_node_map[i].node_id << std::endl;
//                    break;
//            case 764:
//            case 770:
//            case 774:
//            case 778:
//            case 2405555330:
//            case 1418702633:
//            case 289323916:
//                if (used_nodes[i])
//                    std::cout << "+ accessable " << internal_to_external_node_map[i].node_id << ' ' << i << std::endl;
//                else
//                    std::cout << "- unaccessable !!! " << internal_to_external_node_map[i].node_id << std::endl;
//                std::cout << used_nodes[i] << ' ' << std::binary_search(barrier_node_list.begin(), barrier_node_list.end(), i) <<
//                    std::binary_search(traffic_light_list.begin(), traffic_light_list.end(), i) << std::endl;
//                for (unsigned int j: nodes[i]) {
//                    std::cout << j << ' ' << forward_visited[j] << ' ' << backward_visited[j] << std::endl;
//                }
//                for (unsigned int j: restrictions[i]) {
//                    restr = restriction_list[j];
//                    std::cout << restr.fromNode << ' ' << restr.viaNode << ' ' << restr.toNode << ' ' << restr.flags.isOnly << std::endl;
//                }
//        }
//    }

//    for (unsigned int i = 0; i < internal_to_external_node_map.size(); ++i) {
//        if (!used_nodes[i])
//            std::cout << internal_to_external_node_map[i].node_id << ' ' << nodes[i].size() << std::endl;
//        if (internal_to_external_node_map[i].node_id == 359985524) {
//            if (restrictions.find(i) != restrictions.end()) {
//                for (unsigned int j: restrictions[i]) {
//                    std::cout << " from " << internal_to_external_node_map[restriction_list[j].fromNode].node_id <<
//                            " via " << internal_to_external_node_map[restriction_list[j].viaNode].node_id <<
//                            " to " << internal_to_external_node_map[restriction_list[j].toNode].node_id <<
//                            " is_only " << restriction_list[j].flags.isOnly << std::endl;
//                }
//            }
//        }
//    }

    return 0;
}