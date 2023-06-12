#include <iostream>
#include <fstream>
#include "config.h"
#include "solvers/runge_kutta.h"

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

Solver * make_solver(json solver_json, Domain &domain) {
    std::string type = solver_json.at("type");    
    Solver * solver;
    if (type == "runge_kutta") {
	solver = new RungeKutta(solver_json, domain);
    }
    else {
	throw new std::runtime_error("Unknown solver type");
    }
    return solver;
}
