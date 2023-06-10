#include <iostream>
#include <fstream>
#include "config.h"
#include <json.hpp>

using json = nlohmann::json;

Config::Config(std::string file_name) {
    std::ifstream json_file(file_name);
    json json_data = json::parse(json_file);

    this->number_cells = json_data.at("number_cells");
    this->length = json_data.at("length");
    this->velocity = json_data.at("velocity");
    this->max_step = json_data.at("max_step");
    this->max_time = json_data.at("max_time");
    this->print_frequency = json_data.at("print_frequency");
    this->plot_frequency = json_data.at("plot_frequency");
    this->cfl = json_data.at("cfl");
}
