//
// Created by dan1992 on 15.02.16.
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
#include <math.h>

std::string restrictions_path = "map.osrm.restrictions";
std::string path = "map.osrm";

std::string out_restrictions_path = "new_map.osrm.restrictions";
std::string out_path = "new_map.osrm";

std::vector<NodeInfo> internal_to_external_node_map;
std::vector<TurnRestriction> restriction_list;
std::vector<NodeID> barrier_node_list;
std::vector<NodeID> traffic_light_list;
std::vector<ImportEdge> edge_list;

std::vector<bool> visited, unused_nodes;
std::vector<std::set<unsigned int>> nodes;
std::map<NodeID, std::set<TurnRestriction *>> restrictions;

struct Config {
    std::string inputFile;
    std::string outputFile;
};

const int offset = 15;      //point offset, deg*10^6 15

std::ostream& operator<<(std::ostream &s, NodeInfo &n) {
    s << "id " << n.node_id << "\nlat " << n.lat << " lon " << n.lon << std::endl;
    return s;
}

std::ostream& operator<<(std::ostream &s, ImportEdge &e) {
    s << "source " << e.source << " target " << e.target << " name id " << e.name_id << "\nlength " << e.length <<
            " weight " << e.weight << "\nforward " << e.forward << " backward " << e.backward << "\nsplit " <<
            e.is_split << " acces resrticted " << e.access_restricted << " contra flow " << e.contra_flow <<
            " in tiny cc " << e.in_tiny_cc << "\nmaxheight " << e.maxheight << " maxload " << e.maxload << std::endl;
    return s;
}

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
    );

    boost::program_options::variables_map option_variables;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv)
                                          .options(config_options)
                                          .run(),
                                  option_variables);

    boost::program_options::notify(option_variables);

    if (!option_variables.count("inputFile")) {
//        std::cout << "missing input file" << std::endl;
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

NodeID findNextEdge(const NodeID via, const std::vector<bool> &v) {
    if (nodes[via].size() != 2)
        return UINT_MAX;
    for (NodeID n: nodes[via]){
        if (!v[n] && edge_list[n].is_split && edge_list[n].forward && edge_list[n].backward)
            return n;
    }
    return UINT_MAX;
}

std::vector<NodeID> findChain(const NodeID start_edge, std::vector<bool> &v, std::vector<bool> &nv) {
    std::vector<NodeID> chain;
    std::stack<NodeID> stack;
    v[start_edge] = true;
    NodeID via = edge_list[start_edge].source,
            n = findNextEdge(via, v);
    while (n != UINT_MAX) {
        nv[via] = true;
        stack.push(n);
        v[n] = true;
        via = (via == edge_list[n].source) ? edge_list[n].target : edge_list[n].source;
        n = findNextEdge(via, v);
    }
    while (!stack.empty()) {
        chain.push_back(stack.top());
        stack.pop();
    }
    chain.push_back(start_edge);
    via = edge_list[start_edge].target;
    n = findNextEdge(via, v);
    while (n != UINT_MAX) {
        nv[via] = true;
        chain.push_back(n);
        v[n] = true;
        via = (via == edge_list[n].source) ? edge_list[n].target : edge_list[n].source;
        n = findNextEdge(via, v);
    }
//    std::cout << chain.size() << std::endl;
    return chain;
}

float getAngle(const NodeID edge, const NodeID via) {
    int from, to, dx, dy;
    from = via;
    if (edge_list[edge].source == via) {
        to = edge_list[edge].target;
    } else {
        to = edge_list[edge].source;
    }

    NodeInfo to_node = internal_to_external_node_map[to], from_node = internal_to_external_node_map[from];
    ImportEdge *e = &edge_list[edge];

    dx = internal_to_external_node_map[to].lon - internal_to_external_node_map[from].lon;
    dy = internal_to_external_node_map[to].lat - internal_to_external_node_map[from].lat;
    float res;
    if (dx == 0) {
        if (dy == 0)
            return 0;
        res = dy > 0 ? M_PI_2 : -M_PI_2;
    } else {
        res = atan2((float)dy, (float)dx);
    }
    return (res < 0) ? res + M_PI * 2 : res;
}

std::pair<float, float> getClosestAngles(const NodeID edge, const NodeID via) {
    float left = M_PI_2, right = -M_PI_2, curr, temp;
    curr = getAngle(edge, via);
    for (NodeID e: nodes[via]) {
        if (e == edge)
            continue;
        temp = getAngle(e, via) - curr;
        if (temp < -M_PI) {
            temp += M_PI * 2;
        } else
            if (temp > M_PI) {
                temp -= M_PI * 2;
            }
        if (temp > 0) {
            if (left > temp)
                left = temp;
        } else {
            if (right < temp)
                right = temp;
        }
    }
    return std::make_pair(left, right);
}

std::pair<float, float> getBisect(float first, float second) {
    if (second > first) {
        second -= 2 * M_PI;
    }
    std::pair<float, float> res;
    res.first = (first + second) / 2;
    res.second = (res.first > M_PI) ? res.first - M_PI : res.first + M_PI;
    return res;
};

std::pair<int, int> calcOffset(float angle, int d) {
    float lat, lon;
    sincosf(angle, &lat, &lon);
    lon *= d;
    lat *= d;
    return std::make_pair((int)lat, (int)lon);
}

void nodeOffset(const NodeInfo *origin, float angle, int offset, NodeInfo &target) {
    std::pair<int,int> o = calcOffset(angle, offset);
    target.lat = origin->lat + o.first;
    target.lon = origin->lon + o.second;
}

NodeID createNode(const NodeInfo *p, NodeInfo &n, float angle, int offset) {
    nodeOffset(p, angle, offset, n);
    internal_to_external_node_map.push_back(n);
    return internal_to_external_node_map.size() - 1;
}

void createEdge(const ImportEdge *e, NodeID s, NodeID t, bool forward = true, bool sh = true) {
    nodes[s].insert(edge_list.size());
    nodes[t].insert(edge_list.size());
    ImportEdge ie((NodeID)s,
                  (NodeID)t,
                  (NodeID)e->name_id,
                  (EdgeWeight)(sh ? 1 : e->weight),
                  (bool)forward,
                  (bool)!forward,
                  e->type,
                  e->roundabout,
                  e->in_tiny_cc,
                  e->access_restricted,
                  e->contra_flow,
                  (bool)false,
                  (sh ? 1 : e->length),
                  e->maxload,
                  e->maxheight);
    edge_list.push_back(ie);
}

void restrictionReplace(NodeID via, NodeID p_to, NodeID l, NodeID r) {
    if (restrictions.find(via) == restrictions.end())
        return;
    for (TurnRestriction *tr: restrictions[via]) {
        if (tr->fromNode == p_to)
            tr->fromNode = l;
        if (tr->toNode == p_to)
            tr->toNode = r;
    }
}

bool checkTurnRestriction(const NodeID from, const NodeID via, const NodeID to) {
    if (restrictions.find(via) == restrictions.end())
        return true;
    bool allow = true;
    for (TurnRestriction *tr: restrictions[via]) {
        if ((tr->fromNode == from) && (tr->toNode == to))
            return tr->flags.isOnly;
        if (tr->flags.isOnly && (tr->fromNode == from))
            allow = false;
    }
    return allow;
}

void split(std::vector<NodeID> &chain) {
    if (chain.empty())
        return;
    edge_list.reserve(edge_list.size() + (chain.size() + 2 + 1) * 2);
    internal_to_external_node_map.reserve(internal_to_external_node_map.size() + (chain.size() + 1) * 2);
    float angle, s_angle, t_angle, n_angle;
    ImportEdge *edge;
    NodeID edgeID, s_l_i, s_r_i, t_l_i, t_r_i, prev, curr;
    NodeInfo *n_s, *n_t, s_left, s_right, t_left, t_right;
    std::pair<float, float> source, target;
    bool forward = true;

    edgeID = chain[0];
    edge = &edge_list[edgeID];

    if (chain.size() == 1) {
        s_angle = getAngle(edgeID, edge->source);
        t_angle = s_angle > M_PI ? s_angle - M_PI : s_angle + M_PI;

//        std::cout << "one" << std::endl;

//        std::cout << s_angle * 180 / M_PI << " " << t_angle * 180 / M_PI << std::endl;
//        std::cout << *edge << std::endl;

        source = getClosestAngles(edgeID, edge->source);
        target = getClosestAngles(edgeID, edge->target);

//        std::cout << source.first * 180 / M_PI << " " << source.second * 180 / M_PI << std::endl;
//        std::cout << target.first * 180 / M_PI << " " << target.second * 180 / M_PI << std::endl;

        n_s = &internal_to_external_node_map[edge->source];
        n_t = &internal_to_external_node_map[edge->target];

        s_l_i = createNode(n_s, s_left, s_angle + source.first / 2, offset);
        s_r_i = createNode(n_s, s_right, s_angle + source.second / 2, offset);
        t_l_i = createNode(n_t, t_left, t_angle + target.first / 2, offset);
        t_r_i = createNode(n_t, t_right, t_angle + target.second / 2, offset);
        nodes.resize(internal_to_external_node_map.size());

        nodes[edge->source].erase(edgeID);
        nodes[edge->target].erase(edgeID);

        createEdge(edge, edge->source, s_l_i, false);
        createEdge(edge, s_l_i, t_r_i, false, false);
        createEdge(edge, t_r_i, edge->target, false);

        createEdge(edge, edge->source, s_r_i);
        createEdge(edge, s_r_i, t_l_i, true, false);
        createEdge(edge, t_l_i, edge->target);

        if (checkTurnRestriction(edge->target, edge->source, edge->target))
            createEdge(edge, s_l_i, s_r_i);
        if (checkTurnRestriction(edge->source, edge->target, edge->source))
            createEdge(edge, t_l_i, t_r_i);

        restrictionReplace(edge->source, edge->target, s_l_i, s_r_i);
        restrictionReplace(edge->target, edge->source, t_l_i, t_r_i);

//        std::cout << "current edge " << *edge << std::endl;

//        std::cout << "node(s) " << edge->source << " " << internal_to_external_node_map[edge->source] << std::endl;
//        std::cout << "node(l) " << s_l_i << " " << internal_to_external_node_map[s_l_i] << std::endl;
//        std::cout << "node(r) " << s_r_i << " " << internal_to_external_node_map[s_r_i] << std::endl;
//        std::cout << "node(t) " << edge->target << " " << internal_to_external_node_map[edge->target] << std::endl;
//        std::cout << "node(l) " << t_l_i << " " << internal_to_external_node_map[t_l_i] << std::endl;
//        std::cout << "node(r) " << t_r_i << " " << internal_to_external_node_map[t_r_i] << std::endl << std::endl;

//        std::cout << "edge sl " << edge_list[edge_list.size() - 6] << std::endl;
//        std::cout << "edge  l " << edge_list[edge_list.size() - 5] << std::endl;
//        std::cout << "edge tl " << edge_list[edge_list.size() - 4] << std::endl;
//        std::cout << "edge sr " << edge_list[edge_list.size() - 3] << std::endl;
//        std::cout << "edge  r " << edge_list[edge_list.size() - 2] << std::endl;
//        std::cout << "edge tr " << edge_list[edge_list.size() - 1] << std::endl;
        return;
    }

//    std::cout << "chain " << chain.size() << std::endl;

    if (nodes[edge->source].find(chain[1]) != nodes[edge->source].end()) {
        forward = false;
        prev = edge->target;
        curr = edge->source;
    } else {
        prev = edge->source;
        curr = edge->target;
    }
    angle = getAngle(edgeID, prev);
    n_angle = getAngle(chain[1], curr);

//    std::cout << angle * 180 / M_PI << " " << n_angle * 180 / M_PI << std::endl;

    source = getClosestAngles(edgeID, prev);

    n_s = &internal_to_external_node_map[prev];
    n_t = &internal_to_external_node_map[curr];

    s_l_i = createNode(n_s, s_left, angle + source.first / 2, offset);
    s_r_i = createNode(n_s, s_right, angle + source.second / 2, offset);
    nodes.resize(internal_to_external_node_map.size());

    nodes[prev].erase(edgeID);

    angle = angle > M_PI ? angle - M_PI : angle + M_PI;
//    target.first = (angle + n_angle) / 2;
//    target.second = target.first > M_PI ? target.first - M_PI : target.first + M_PI;

    target = getBisect(angle, n_angle);

    t_l_i = createNode(n_t, t_left, target.second, offset);
    t_r_i = createNode(n_t, t_right, target.first, offset);
    nodes.resize(internal_to_external_node_map.size());

    createEdge(edge, prev, s_l_i, false);
    createEdge(edge, s_l_i, t_r_i, false, false);

    createEdge(edge, prev, s_r_i);
    createEdge(edge, s_r_i, t_l_i, true, false);

    if (checkTurnRestriction(curr, prev, curr))
        createEdge(edge, s_l_i, s_r_i);

    restrictionReplace(prev, curr, s_l_i, s_r_i);

//    std::cout << "\n\ncurrent edge " << *edge << std::endl;

//    std::cout << angle * 180 / M_PI << " " << n_angle * 180 / M_PI << std::endl;

//    std::cout << source.first * 180 / M_PI << " " << source.second * 180 / M_PI << std::endl;
//    std::cout << target.first * 180 / M_PI << " " << target.second * 180 / M_PI << std::endl;

//    std::cout << "node(s) " << prev << " " << internal_to_external_node_map[prev] << std::endl;
//    std::cout << "node(l) " << s_l_i << " " << internal_to_external_node_map[s_l_i] << std::endl;
//    std::cout << "node(r) " << s_r_i << " " << internal_to_external_node_map[s_r_i] << std::endl;
//    std::cout << "node(t) " << curr << " " << internal_to_external_node_map[curr] << std::endl;
//    std::cout << "node(l) " << t_l_i << " " << internal_to_external_node_map[t_l_i] << std::endl;
//    std::cout << "node(r) " << t_r_i << " " << internal_to_external_node_map[t_r_i] << std::endl << std::endl;

//    std::cout << "edge sl " << edge_list[edge_list.size() - 4] << std::endl;
//    std::cout << "edge  l " << edge_list[edge_list.size() - 3] << std::endl;
//    std::cout << "edge sr " << edge_list[edge_list.size() - 2] << std::endl;
//    std::cout << "edge  r " << edge_list[edge_list.size() - 1] << std::endl;

    for (int i = 1; i < chain.size() - 1; i++){
        s_l_i = t_r_i;
        s_r_i = t_l_i;

        edgeID = chain[i];
        edge = &edge_list[edgeID];

        prev = curr;
        curr = (prev == edge->source) ? edge->target : edge->source;

        n_t = &internal_to_external_node_map[curr];

        angle = n_angle > M_PI ? n_angle - M_PI : n_angle + M_PI;
        n_angle = getAngle(chain[i + 1], curr);
//        target.first = (angle + n_angle) / 2;
//        target.second = target.first > M_PI ? target.first - M_PI : target.first + M_PI;

        target = getBisect(angle, n_angle);

        t_l_i = createNode(n_t, t_left, target.second, offset);
        t_r_i = createNode(n_t, t_right, target.first, offset);
        nodes.resize(internal_to_external_node_map.size());

        createEdge(edge, s_l_i, t_r_i, false, false);

        createEdge(edge, s_r_i, t_l_i, true, false);

//        std::cout << "current edge " << *edge << std::endl;

//        std::cout << angle * 180 / M_PI << " " << n_angle * 180 / M_PI << std::endl;

//        std::cout << source.first * 180 / M_PI << " " << source.second * 180 / M_PI << std::endl;
//        std::cout << target.first * 180 / M_PI << " " << target.second * 180 / M_PI << std::endl;

//        std::cout << "node(t) " << curr << " " << internal_to_external_node_map[curr] << std::endl;
//        std::cout << "node(l) " << t_l_i << " " << internal_to_external_node_map[t_l_i] << std::endl;
//        std::cout << "node(r) " << t_r_i << " " << internal_to_external_node_map[t_r_i] << std::endl << std::endl;

//        std::cout << "edge  l " << edge_list[edge_list.size() - 2] << std::endl;
//        std::cout << "edge  r " << edge_list[edge_list.size() - 1] << std::endl;
    }

    s_l_i = t_r_i;
    s_r_i = t_l_i;

    edgeID = chain.back();
    edge = &edge_list[edgeID];

    prev = curr;
    curr = (prev == edge->source) ? edge->target : edge->source;

    n_t = &internal_to_external_node_map[curr];

    n_angle = n_angle > M_PI ? n_angle - M_PI : n_angle + M_PI;
    target = getClosestAngles(chain.back(), curr);

    t_l_i = createNode(n_t, t_left, n_angle + target.first / 2, offset);
    t_r_i = createNode(n_t, t_right, n_angle + target.second / 2, offset);
    nodes.resize(internal_to_external_node_map.size());

    createEdge(edge, s_l_i, t_r_i, false, false);
    createEdge(edge, t_r_i, curr, false);

    createEdge(edge, s_r_i, t_l_i, true, false);
    createEdge(edge, t_l_i, curr);

    if (checkTurnRestriction(prev, curr, prev))
        createEdge(edge, t_l_i, t_r_i);


    restrictionReplace(curr, prev, t_l_i, t_r_i);

//    std::cout << "current edge " << *edge << std::endl;

//    std::cout << angle * 180 / M_PI << " " << n_angle * 180 / M_PI << std::endl;

//    std::cout << source.first * 180 / M_PI << " " << source.second * 180 / M_PI << std::endl;
//    std::cout << target.first * 180 / M_PI << " " << target.second * 180 / M_PI << std::endl;

//    std::cout << "node(s) " << edge->source << " " << internal_to_external_node_map[edge->source] << std::endl;
//    std::cout << "node(l) " << s_l_i << " " << internal_to_external_node_map[s_l_i] << std::endl;
//    std::cout << "node(r) " << s_r_i << " " << internal_to_external_node_map[s_r_i] << std::endl;
//    std::cout << "node(t) " << edge->target << " " << internal_to_external_node_map[edge->target] << std::endl;
//    std::cout << "node(l) " << t_l_i << " " << internal_to_external_node_map[t_l_i] << std::endl;
//    std::cout << "node(r) " << t_r_i << " " << internal_to_external_node_map[t_r_i] << std::endl << std::endl;

//    std::cout << "edge  l " << edge_list[edge_list.size() - 4] << std::endl;
//    std::cout << "edge tl " << edge_list[edge_list.size() - 3] << std::endl;
//    std::cout << "edge  r " << edge_list[edge_list.size() - 2] << std::endl;
//    std::cout << "edge tr " << edge_list[edge_list.size() - 1] << std::endl;
    return;
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
        restrictions[restriction_list[i].viaNode].insert(&restriction_list[i]);
    }

//    unsigned s = 0, c = 0, z = 0;
//    for (ImportEdge &e: edge_list) {
//        NodeInfo sn = internal_to_external_node_map[e.source], tn = internal_to_external_node_map[e.target];
//        if (sn.lat == tn.lat && sn.lon == tn.lon) {
//            if (sn.lat < 50610000 && sn.lat > 50300000 && sn.lon < 30710000 && sn.lon > 30290000) {
//                std::cout << sn.node_id << " = " << tn.node_id << std::endl;
//                z++;
//            }
//        }
//    }
//    std::cout << "coliding nodes " << 2 * z << std::endl;

    visited.resize(edge_list.size(), false);
    unused_nodes.resize(internal_to_external_node_map.size(), false);
    std::vector<std::vector<NodeID>> chains;
    unsigned long edges_s = edge_list.size();
    ImportEdge *edge;
    NodeID t = 0;

//    for (ImportEdge &e: edge_list) {
//        if (internal_to_external_node_map[e.source].node_id == 1745493511 ||
//            internal_to_external_node_map[e.target].node_id == 1745493511) {
//            t = e.name_id;
//            break;
//        }
//    }
//    for (ImportEdge &e: edge_list) {
//        if (e.name_id == t) {
//            std::cout << "street " <<  internal_to_external_node_map[e.source].node_id << " " <<
//                    internal_to_external_node_map[e.target].node_id << std::endl;
//            e.is_split = true;
//        }
//    }


    for (NodeID i = 0; i < edges_s; i++) {
        edge = &edge_list[i];
        if (visited[i] || !edge_list[i].is_split || !edge_list[i].forward || !edge_list[i].backward)
            continue;
        chains.push_back(findChain(i, visited, unused_nodes));
    }

    for (std::vector<NodeID> &chain: chains) {
        split(chain);
    }

    visited.resize(edge_list.size(), false);
    unused_nodes.resize(internal_to_external_node_map.size(), false);

    NodeID new_nodes = 0, new_edges = 0;
    for (unsigned int i = 0; i < n; ++i) {
        if (!unused_nodes[i]) new_nodes++;
    }
    std::cout << n << ' ' << new_nodes << std::endl;

    n = UINT_MAX - 1;
    for (int i = internal_to_external_node_map.size() - 1; i >= 0; i--) {
        if (internal_to_external_node_map[i].node_id == UINT_MAX) {
            internal_to_external_node_map[i].node_id = n--;
        }
    }

    std::ofstream out_restrictions(out_restrictions_path, std::ios::binary);
    out_restrictions.write((char *) &fingerPrint, sizeof(FingerPrint));
    out_restrictions.write((char *) &n, sizeof(unsigned));

    unsigned int used_restrictions = 0;

    for (TurnRestriction restriction : restriction_list) {
        if (!unused_nodes[restriction.viaNode] &&
            !unused_nodes[restriction.fromNode] &&
            !unused_nodes[restriction.toNode]) {
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
        if (!unused_nodes[i]) {
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
        if (!visited[i]) {
            new_edges++;
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
    std::ios::pos_type current = out.tellp();
    out.seekp(previous);
    out.write((char *) &new_edges, sizeof(unsigned));
    out.seekp(current);

    std::cout << new_edges << std::endl;

    out.close();

//    visited.clear();

//    std::pair<float, float> pair;
//    std::pair<int, int> o;
//    float a, min = M_PI_2;
//
//    for (int i = 0; i < 100; i++) {
//        std::cout << "\n\nnode " << internal_to_external_node_map[i].node_id << std::endl;
//        for (NodeID e: nodes[i]) {
//            a = getAngle(e, i);
//            pair = getClosestAngles(e, i);
//            if (pair.first > 0 && pair.first < min)
//                min = pair.first;
//            if (pair.second > 0 && pair.second < min)
//                min = pair.second;
//            o = calcOffset(a, offset);
//            std::cout << "\nangle " << a * 180 / M_PI << "\nleft " << pair.first * 180 / M_PI << "\nright " <<
//                    pair.second * 180 / M_PI <<  "\nlat " << o.first << "\nlon " << o.second << std::endl;
//        }
//    }

//    std::cout << s << std::endl << c << std::endl;
    return 0;
}