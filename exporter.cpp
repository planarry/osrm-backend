

#include "Util/GitDescription.h"
#include "Util/ProgramOptions.h"
#include "Util/SimpleLogger.h"
#include "Util/TimingUtil.h"
#include "Util/FingerPrint.h"
#include "DataStructures/Restriction.h"
#include "DataStructures/ImportNode.h"
#include "DataStructures/RangeTable.h"
#include "DataStructures/Percent.h"
#include "typedefs.h"

#include <vector>
#include <map>
#include <chrono>
#include <fstream>
#include <iostream>
#include <pqxx/pqxx>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

struct Config {
  std::string fileName;
  std::string dbHost;
  std::string dbUser;
  std::string dbName;
  std::string dbPass;
  boost::filesystem::path ini;
};

struct Edge {
  unsigned source;
  unsigned target;
  int length;
  short dir;
  int weight;
  short type;
  unsigned nameID;
  bool is_roundabout;
  bool ignore_in_grid;
  bool is_access_restricted;
  bool is_contra_flow;
  bool is_split;
  short maxload;
  short maxheight;
};
const unsigned int PACK_SIZE=1000;

bool ParseArguments(int argc, char *argv[], Config& config)
{
  // declare a group of options that will be allowed only on command line
  boost::program_options::options_description generic_options("Options");
  generic_options.add_options()("version,v", "Show version")("help,h", "Show this help message")(
      "config,c",
      boost::program_options::value<boost::filesystem::path>(&config.ini)
        ->default_value("exporter.ini"),
      "Path to a configuration file.");

  // declare a group of options that will be allowed both on command line and in config file
  boost::program_options::options_description config_options("Configuration");
  config_options.add_options()(
        "fileName,o",
        boost::program_options::value<std::string>(&config.fileName),
        "Outpup *.osrm filename"
      )(
        "dbHost,H",
        boost::program_options::value<std::string>(&config.dbHost)->default_value("localhost"),
        "PostrgeSQL host"
      )(
        "dbUser,u",
        boost::program_options::value<std::string>(&config.dbUser)->default_value("postgres"),
        "PostrgeSQL user"
      )(
        "dbPass,p",
        boost::program_options::value<std::string>(&config.dbPass)->default_value("postgres"),
        "PostrgeSQL password"
      )(
        "dbName,b",
        boost::program_options::value<std::string>(&config.dbName),
        "PostrgeSQL input database name"
      );

  // combine above options for parsing
  boost::program_options::options_description cmdline_options;
  cmdline_options.add(generic_options).add(config_options);

  boost::program_options::options_description config_file_options;
  config_file_options.add(config_options);

  boost::program_options::options_description visible_options(
      boost::filesystem::basename(argv[0]) + " [options]");
  visible_options.add(generic_options).add(config_options);

  // parse command line options
  boost::program_options::variables_map option_variables;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv)
                                    .options(cmdline_options)
                                    .run(),
                                option_variables);

  if (option_variables.count("version"))
  {
      SimpleLogger().Write() << g_GIT_DESCRIPTION;
      return false;
  }

  if (option_variables.count("help"))
  {
      SimpleLogger().Write() << visible_options;
      return false;
  }

  boost::program_options::notify(option_variables);

  // parse config file
  if (boost::filesystem::is_regular_file(config.ini))
  {
      SimpleLogger().Write() << "Reading options from: " << config.ini.string();
      std::string ini_file_contents = ReadIniFileAndLowerContents(config.ini);
      std::stringstream config_stream(ini_file_contents);
      boost::program_options::store(parse_config_file(config_stream, config_file_options),
                                    option_variables);
      boost::program_options::notify(option_variables);
  }

  if (!option_variables.count("dbName") || !option_variables.count("fileName"))
  {
      SimpleLogger().Write() << visible_options;
      return false;
  }

  return true;
}

int main (int argc, char *argv[])
{
  Config config;
  const FingerPrint fingerprint;
  std::vector<std::string> names_list;
  std::map<std::string, unsigned int> string_map;
  
  LogPolicy::GetInstance().Unmute();
  if(!ParseArguments(argc, argv, config))
    return 0;
  
  std::string conStr("");
  conStr.append(" dbname="   + config.dbName);
  conStr.append(" host="     + config.dbHost);
  conStr.append(" user="     + config.dbUser);
  conStr.append(" password=" + config.dbPass);

  std::shared_ptr<pqxx::connection> con;
  try {
      con = std::make_shared<pqxx::connection>(conStr);
  } catch (const pqxx::sql_error &e) {
      SimpleLogger().Write(logWARNING) << e.what();
      return 1;
  }
  pqxx::work w(*con);
  pqxx::result res;
        
  unsigned nodes_count;
  unsigned edges_count;
  unsigned turns_count;
  SimpleLogger().Write() << "Fetching counts...";
  TIMER_START(counts);
  try {
    nodes_count=w.exec("SELECT count(*) FROM nodes WHERE confirmed").front()[0].as< unsigned int >();
  } catch (const std::exception &e) {
    SimpleLogger().Write(logWARNING) << e.what();
  }   
  try {
    edges_count=w.exec("SELECT count(*) FROM edges e "
      "INNER JOIN nodes s on srcID=s.ID "
      "INNER JOIN nodes t on trgID=t.ID "
      "WHERE e.confirmed AND s.confirmed AND t.confirmed "
      "AND speed>0 AND (forward or backward) "
      "AND ST_Distance(s.coord::geometry, t.coord::geometry)>0").front()[0].as< unsigned int >();
  } catch (const std::exception &e) {
    SimpleLogger().Write(logWARNING) << e.what();
  } 
  try {
    turns_count=w.exec("SELECT count(*) FROM turns r "
      "INNER JOIN nodes s on srcID=s.ID "
      "INNER JOIN nodes v on viaID=v.ID "
      "INNER JOIN nodes t on trgID=t.ID "
      "WHERE r.confirmed AND s.confirmed AND v.confirmed AND t.confirmed").front()[0].as< unsigned int >();
  } catch (const std::exception &e) {
    SimpleLogger().Write(logWARNING) << e.what();
  } 
  TIMER_STOP(counts);
  SimpleLogger().Write() << "ok after " << TIMER_SEC(counts) << "s";
  SimpleLogger().Write() << "Nodes: " << nodes_count;
  SimpleLogger().Write() << "Edges: " << edges_count;
  SimpleLogger().Write() << "Turns: " << turns_count;
  Percent progress(nodes_count+edges_count+turns_count);
  
  
  SimpleLogger().Write() << "Fetching and writing files...";
  TIMER_START(write);
  std::ofstream file_out_stream;
  file_out_stream.open(config.fileName.c_str(), std::ios::binary);
  file_out_stream.write((char *)&fingerprint, sizeof(FingerPrint));
  file_out_stream.write((char *)&nodes_count, sizeof(unsigned));
  try {
    auto query = "SELECT ID, "
      "(ST_Y(coord::geometry)*1e6)::int as lat, "
      "(ST_X(coord::geometry)*1e6)::int as lon, "
      "trafficlight "
      "FROM nodes "
      "WHERE confirmed "
      "ORDER BY ID";
    pqxx::icursorstream cur(w, query, "cur", PACK_SIZE);
    while(cur>>res)
      for(auto row : res)
      {
        progress.printIncrement();
        ExternalMemoryNode node;
        node.lat = row["lat"].as< int >();
        node.lon = row["lon"].as< int >();
        node.node_id = row["ID"].as< unsigned int >();
        node.bollard = 0;
        node.trafficLight = row["trafficlight"].as< bool >();
        file_out_stream.write((char *)&(node), sizeof(ExternalMemoryNode));
      }
  } catch (const std::exception &e) {
    SimpleLogger().Write(logWARNING) << e.what();
  }   
  
  
  file_out_stream.write((char *)&edges_count, sizeof(unsigned));
  try {
    auto query = "SELECT srcID, trgID, "
      "forward, backward, "
      "splitted, streetname, speed, "
      //"ST_Distance(s.coord, t.coord)::int length, "
      "(ST_Y(s.coord::geometry)*1e6)::int as lat1, "
      "(ST_X(s.coord::geometry)*1e6)::int as lon1, "
      "(ST_Y(t.coord::geometry)*1e6)::int as lat2, "
      "(ST_X(t.coord::geometry)*1e6)::int as lon2, "
      "(ST_Distance(s.coord, t.coord)/(speed/3.6)*10)::int weight, "
      "coalesce(nullif(maxload,0)/100, 255)::int maxload, " //max 25500 kilograms
      "coalesce(nullif(maxheight,0)*10, 63)::int maxheight " //max 6.3 meters height
      "FROM edges e "
      "INNER JOIN nodes s on srcID=s.ID "
      "INNER JOIN nodes t on trgID=t.ID "
      "WHERE e.confirmed AND s.confirmed AND t.confirmed "
      "AND speed>0 "
      "AND (forward or backward) "
      "AND ST_Distance(s.coord, t.coord)>0 "
      "ORDER BY srcID";
    pqxx::icursorstream cur(w, query, "cur", PACK_SIZE);
    while(cur>>res)
      for(auto row : res)
      {
        progress.printIncrement();
        Edge edge;
        edge.source = row["srcID"].as< unsigned >();
        edge.target = row["trgID"].as< unsigned >();
        edge.length = FixedPointCoordinate::ApproximateDistance(
          row["lat1"].as<int>(), row["lon1"].as<int>(),
          row["lat2"].as<int>(), row["lon2"].as<int>());
        //edge.length = std::max(1, row["length"].as< int >());
        edge.dir = row["forward"].as< bool >() ? row["forward"].as< bool >()!=row["backward"].as< bool >() : 2;
        edge.weight = std::max(1, row["weight"].as< int >());
        edge.type = 0;
        edge.is_roundabout = 0;
        edge.ignore_in_grid = 0;
        edge.is_access_restricted = 0;
        edge.is_contra_flow = 0;
        edge.is_split = row["splitted"].as< bool >();
        edge.maxload = row["maxload"].as< short >();
        edge.maxheight = row["maxheight"].as< short >();
        
        // Get the unique identifier for the street name
        auto name = row["streetname"].as< std::string >();
        if(!name.empty())
        {
          const auto &string_map_iterator = string_map.find(name);
          if (string_map.end() == string_map_iterator)
          {
            edge.nameID = names_list.size();
            names_list.push_back(name);
            string_map.insert(std::make_pair(name, edge.nameID));
          }
          else
            edge.nameID = string_map_iterator->second;
        }
        else edge.nameID = 0;
        
        file_out_stream.write((char *)&edge.source, sizeof(unsigned));
        file_out_stream.write((char *)&edge.target, sizeof(unsigned));
        file_out_stream.write((char *)&edge.length, sizeof(int));
        file_out_stream.write((char *)&edge.dir, sizeof(short));
        file_out_stream.write((char *)&edge.weight, sizeof(int));
        file_out_stream.write((char *)&edge.type, sizeof(short));
        file_out_stream.write((char *)&edge.nameID, sizeof(unsigned));
        file_out_stream.write((char *)&edge.is_roundabout, sizeof(bool));
        file_out_stream.write((char *)&edge.ignore_in_grid, sizeof(bool));
        file_out_stream.write((char *)&edge.is_access_restricted, sizeof(bool));
        file_out_stream.write((char *)&edge.is_contra_flow, sizeof(bool));
        file_out_stream.write((char *)&edge.is_split, sizeof(bool));
        file_out_stream.write((char *)&edge.maxload, sizeof(short));
        file_out_stream.write((char *)&edge.maxheight, sizeof(short));
      }
  } catch (const std::exception &e) {
    SimpleLogger().Write(logWARNING) << e.what();
  }  
  file_out_stream.close();
  
  
  std::ofstream restrictions_out_stream;
  restrictions_out_stream.open((config.fileName+".restrictions").c_str(), std::ios::binary);
  restrictions_out_stream.write((char *)&fingerprint, sizeof(FingerPrint));
  restrictions_out_stream.write((char *)&turns_count, sizeof(unsigned));
  try {
    auto query = "SELECT srcID, viaID, trgID "
      "FROM turns r "
      "INNER JOIN nodes s on srcID=s.ID "
      "INNER JOIN nodes v on viaID=v.ID "
      "INNER JOIN nodes t on trgID=t.ID "
      "WHERE r.confirmed AND s.confirmed AND v.confirmed AND t.confirmed "
      "ORDER BY srcID";
    pqxx::icursorstream cur(w, query, "cur", PACK_SIZE);
    while(cur>>res)
      for(auto row : res)
      {
        progress.printIncrement();
        TurnRestriction turn;
        turn.fromNode = row["srcID"].as< unsigned int >();
        turn.viaNode = row["viaID"].as< unsigned int >();
        turn.toNode = row["trgID"].as< unsigned int >();
        restrictions_out_stream.write((char *)&(turn), sizeof(TurnRestriction));
      }
  } catch (const std::exception &e) {
    SimpleLogger().Write(logWARNING) << e.what();
  }  
  restrictions_out_stream.close();
  
  boost::filesystem::ofstream name_file_stream(config.fileName+".names", std::ios::binary);
  unsigned total_length = 0;
  std::vector<unsigned> name_lengths;
  for (const std::string &temp_string : names_list)
  {
      const unsigned string_length = std::min(static_cast<unsigned>(temp_string.length()), 255u);
      name_lengths.push_back(string_length);
      total_length += string_length;
  }
  RangeTable<> table(name_lengths);
  name_file_stream << table;
  name_file_stream.write((char*) &total_length, sizeof(unsigned));
  for (const std::string &temp_string : names_list)
  {
      const unsigned string_length = std::min(static_cast<unsigned>(temp_string.length()), 255u);
      name_file_stream.write(temp_string.c_str(), string_length);
  }
  name_file_stream.close();
  
  TIMER_STOP(write);
  SimpleLogger().Write() << "ok after " << TIMER_SEC(write) << "s";
}
