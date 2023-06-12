#ifndef _CONFIG_H
#define _CONFIG_H

#include <string>
#include <json.hpp>
#include "solvers/solver.h"
#include "domain.h"

using json = nlohmann::json;

struct Config {
public:
    Config(std::string);
    json domain_json();
    json equation_json();
    json solver_json();

private:
    std::string _title;
    json _json_data;
};

Solver * make_solver(json solver_json, Domain &domain);

#endif
