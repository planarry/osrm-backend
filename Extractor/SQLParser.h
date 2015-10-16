//
// Created by samogot on 07.10.15.
//

#ifndef OSRM_SQLPARSER_H
#define OSRM_SQLPARSER_H


#include "BaseParser.h"
#include "../DataStructures/Restriction.h"
#include <pqxx/pqxx>
#include <memory>


class ExtractorCallbacks;

class SQLParser : public BaseParser {

    std::shared_ptr <pqxx::connection> connection;
    const char *connection_string;

public:
    SQLParser(const char *connection_string,
              ExtractorCallbacks *extractor_callbacks,
              ScriptingEnvironment &scripting_environment);

    bool ReadHeader();

    bool Parse();

};


#endif //OSRM_SQLPARSER_H
