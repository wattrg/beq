#include <iostream>
#include <fstream>
#include "config.h"

Config::Config(std::string file_name) {
    std::ifstream json_file(file_name);
    json json_data = json::parse(json_file);

    _title = json_data.at("title");
    _json_data = json_data;
}

json Config::domain_json() {
    return _json_data.at("domain");
}

json Config::equation_json() {
    return _json_data.at("equation");
}

json Config::solver_json() {
    return _json_data.at("solver");
}

json Config::gas_model() {
    return _json_data.at("gas_model");
}

