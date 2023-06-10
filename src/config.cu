#include <iostream>
#include <fstream>
#include "config.h"
#include <json.hpp>

using json = nlohmann::json;

Config::Config(std::string file_name) {
    std::ifstream json_file(file_name);
    json json_data = json::parse(json_file);

    this->_number_cells = json_data.at("number_cells");
    this->_length = json_data.at("length");
    this->_velocity = json_data.at("velocity");
    this->_max_step = json_data.at("max_step");
    this->_max_time = json_data.at("max_time");
    this->_print_frequency = json_data.at("print_frequency");
    this->_plot_frequency = json_data.at("plot_frequency");
    this->_cfl = json_data.at("cfl");
}
