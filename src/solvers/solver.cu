#include "solver.h"
#include "runge_kutta.h"

Solver * make_solver(json solver_json, Domain &domain, Equation *equation) {
    std::string type = solver_json.at("type");    
    Solver * solver;
    if (type == "runge_kutta") {
	solver = new RungeKutta(solver_json, domain, equation);
    }
    else {
	throw new std::runtime_error("Unknown solver type");
    }
    return solver;
}
