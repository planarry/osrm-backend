//
// Created by samogot on 07.10.15.
//

#include "SQLParser.h"

#include "ExtractionWay.h"
#include "ExtractorCallbacks.h"

#include "../DataStructures/HashTable.h"
#include "../DataStructures/ImportNode.h"
#include "../DataStructures/Restriction.h"
#include "../Util/SimpleLogger.h"
#include "../DataStructures/Percent.h"
#include "../Util/StringUtil.h"

SQLParser::SQLParser(const char *connection_string,
                     ExtractorCallbacks *extractor_callbacks,
                     ScriptingEnvironment &scripting_environment)
        : BaseParser(extractor_callbacks, scripting_environment), connection_string(connection_string) { }

bool SQLParser::ReadHeader() {
    try {
        connection = std::make_shared<pqxx::connection>(connection_string);
        connection->set_variable("application_name", "'osrm-extractor'");
        return true;
    } catch (const pqxx::failure &e) {
        SimpleLogger().Write(logWARNING) << e.what();
        return false;
    }
}

const unsigned int PACK_SIZE = 10000;

bool SQLParser::Parse() {
    pqxx::work w(*connection);
    pqxx::result res, inner_res;
    unsigned nodes_count, ways_count, relations_count;
    Percent progress(0);

    try {
        nodes_count = (unsigned) w.exec(
                "SELECT reltuples FROM pg_class WHERE relname = 'current_nodes'").front()[0].as<double>();
        ways_count = (unsigned) w.exec(
                "SELECT reltuples FROM pg_class WHERE relname = 'current_ways'").front()[0].as<double>();
        relations_count = (unsigned) w.exec("SELECT count(distinct relation_id) "
                                                    "FROM current_relation_members "
                                                    "WHERE member_role='from' AND member_type='Way'").front()[0].as<double>();
    } catch (const std::exception &e) {
        SimpleLogger().Write(logWARNING) << e.what();
        return false;
    }


    SimpleLogger().Write(logINFO) << "Fetching ~" << nodes_count << " nodes:";
    progress.reinit(nodes_count, 1);
    try {
        auto query = "SELECT id, round(latitude/10)::int lat, round(longitude/10)::int lon, v, k "
                "FROM current_nodes n "
                "LEFT JOIN current_node_tags ON node_id=n.id "
                "ORDER BY id";
        pqxx::icursorstream cur(w, query, "cur", PACK_SIZE);
        unsigned prev_node_id = UINT_MAX;
        ImportNode node;

        while (cur >> res)
            for (const auto &row : res) {
                if (row["id"].as<unsigned>() != prev_node_id) {
                    if (prev_node_id != UINT_MAX) {
                        ParseNodeInLua(node, lua_state);
                        extractor_callbacks->ProcessNode(node);
                        progress.printIncrement();
                        node = ImportNode();
                    }
                    node.lat = row["lat"].as<int>();
                    node.lon = row["lon"].as<int>();
                    prev_node_id = node.node_id = row["id"].as<unsigned>();
                    if (!row["k"].is_null())
                        node.keyVals.Add(row["k"].as<std::string>(), row["v"].as<std::string>());
                }
                else
                    node.keyVals.Add(row["k"].as<std::string>(), row["v"].as<std::string>());
            }
    } catch (const std::exception &e) {
        SimpleLogger().Write(logWARNING) << e.what();
        return false;
    }


    SimpleLogger().Write(logINFO) << "Fetching ~" << ways_count << " ways:";
    progress.reinit(ways_count, 1);
    try {
        auto query = "SELECT id, k, v,"
                "trim(BOTH '{}' FROM (SELECT array_agg(node_id ORDER BY sequence_id) FROM current_way_nodes WHERE way_id=w.id)::varchar) nodes,"
                "(SELECT avg((SELECT v FROM current_relation_tags rt WHERE rmf.relation_id=rt.relation_id AND k='speed')::int)"
                " FROM current_way_nodes n1"
                " INNER JOIN current_way_nodes n2 ON (n1.way_id=n2.way_id AND n1.sequence_id+1=n2.sequence_id)"
                " INNER JOIN current_relation_members rmf ON (rmf.member_type='Node' AND rmf.member_role='from' AND rmf.member_id=n2.node_id)"
                " INNER JOIN current_relation_members rmt ON (rmt.member_type='Node' AND rmt.member_role='to' AND rmt.member_id=n1.node_id AND rmf.relation_id=rmt.relation_id)"
                " WHERE n1.way_id=w.id) as b_speed,"
                "(SELECT avg((select v from current_relation_tags rt where rmf.relation_id=rt.relation_id and k='speed')::int)"
                " FROM current_way_nodes n1"
                " INNER JOIN current_way_nodes n2 ON (n1.way_id=n2.way_id AND n1.sequence_id+1=n2.sequence_id)"
                " INNER JOIN current_relation_members rmf ON (rmf.member_type='Node' AND rmf.member_role='from' AND rmf.member_id=n1.node_id)"
                " INNER JOIN current_relation_members rmt ON (rmt.member_type='Node' AND rmt.member_role='to' AND rmt.member_id=n2.node_id AND rmf.relation_id=rmt.relation_id)"
                " WHERE n1.way_id=w.id) as f_speed "
                "FROM current_ways w "
                "LEFT JOIN current_way_tags ON way_id=w.id ";
        pqxx::icursorstream cur(w, query, "cur", PACK_SIZE);
        unsigned prev_way_id = UINT_MAX;
        ExtractionWay way;
        double temp_speed = -1, temp_backward_speed = -1;

        while (cur >> res)
            for (auto row : res) {
                if (row["id"].as<unsigned>() != prev_way_id) {
                    if (prev_way_id != UINT_MAX) {
                        ParseWayInLua(way, lua_state);
                        if (way.speed > 0) {
                            if (temp_speed > 0)
                                way.speed = temp_speed;
                            if (temp_backward_speed > 0)
                                way.backward_speed = temp_backward_speed;
                        }
                        extractor_callbacks->ProcessWay(way);
                        progress.printIncrement();
                        way = ExtractionWay();
                        temp_speed = -1;
                        temp_backward_speed = -1;
                    }
                    prev_way_id = way.id = row["id"].as<unsigned>();
                    if (!row["b_speed"].is_null()) {
                        temp_speed = row["b_speed"].as<double>();
                        temp_backward_speed = row["b_speed"].as<double>();
                    }
                    if (!row["f_speed"].is_null())
                        temp_speed = row["f_speed"].as<double>();
                    if (!row["nodes"].is_null()) {
                        std::stringstream ss(row["nodes"].as<std::string>());
                        std::string node_id_str;
                        while (std::getline(ss, node_id_str, ','))
                            way.path.push_back(StringToUint(node_id_str));
                    }
                    if (!row["k"].is_null())
                        way.keyVals.Add(row["k"].as<std::string>(), row["v"].as<std::string>());
                }
                else
                    way.keyVals.Add(row["k"].as<std::string>(), row["v"].as<std::string>());
            }
    } catch (const std::exception &e) {
        SimpleLogger().Write(logWARNING) << e.what();
        return false;
    }


    if (use_turn_restrictions) {
        SimpleLogger().Write(logINFO) << "Fetching " << relations_count << " relations:";
        progress.reinit(relations_count, 5);
        try {
            auto query = "SELECT "
                    "(SELECT v FROM current_relation_tags"
                    " WHERE relation_id = r.id"
                    " AND k='except' LIMIT 1) as except,"
                    "(SELECT 1 FROM current_relation_tags"
                    " WHERE relation_id = r.id"
                    " AND k='restriction' AND v LIKE 'only\\_%' LIMIT 1) IS NOT NULL AS is_only,"
                    "(SELECT member_id FROM current_relation_members"
                    " WHERE relation_id = r.id"
                    " AND member_role='from' AND member_type='Way' LIMIT 1) as from,"
                    "(SELECT member_id FROM current_relation_members"
                    " WHERE relation_id = r.id"
                    " AND member_role='to' AND member_type='Way' LIMIT 1) as to,"
                    "(SELECT member_id FROM current_relation_members"
                    " WHERE relation_id = r.id"
                    " AND member_role='via' AND member_type='Node' LIMIT 1) as via "
                    "FROM current_relations r "
                    "WHERE (SELECT 1 FROM current_relation_members"
                    " WHERE relation_id = r.id"
                    " AND member_role='from' AND member_type='Way' LIMIT 1) IS NOT NULL"
                    " AND (SELECT 1 FROM current_relation_tags"
                    " WHERE relation_id = r.id"
                    " AND k='type' AND v='restriction' LIMIT 1) IS NOT NULL";
            pqxx::icursorstream cur(w, query, "cur", PACK_SIZE);
            while (cur >> res)
                for (const auto &row : res) {
                    progress.printIncrement();
                    if (!row["except"].is_null()
                        && ShouldIgnoreRestriction(row["except"].as<std::string>()))
                        continue;
                    InputRestrictionContainer restriction;
                    restriction.restriction.flags.isOnly = row["is_only"].as<bool>();
                    if (!row["from"].is_null())
                        restriction.fromWay = row["from"].as<unsigned>();
                    if (!row["to"].is_null())
                        restriction.toWay = row["to"].as<unsigned>();
                    if (!row["via"].is_null())
                        restriction.restriction.viaNode = row["via"].as<unsigned>();
                    extractor_callbacks->ProcessRestriction(restriction);
                }
        } catch (const std::exception &e) {
            SimpleLogger().Write(logWARNING) << e.what();
            return false;
        }
    }
    return true;
}
