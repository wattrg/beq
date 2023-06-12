#include <exception>
#include <iostream>
#include <fstream>
#include <math.h>
#include <filesystem>
#include <stdexcept>

#include "config.h"
#include "domain.h"
#include "equations/equation.h"
#include "solvers/solver.h"
#include "solvers/runge_kutta.h"


#define _STRINGIFY(x) #x
#define STRINGIFY(x) _STRINGIFY(x)

void print_header() {
    std::cout << "beq: Boltzmann equation solver\n";
    std::cout << "Git branch: " << STRINGIFY(GIT_BRANCH) << "\n";
    std::cout << "Git commit: " << STRINGIFY(GIT_HASH) << "\n";
    std::cout << "Build date: " << STRINGIFY(COMPILE_TIME) << "\n";
}


void read_initial_condition(double *phi, int n) {
    std::ifstream initial_condition("solution/phi_0.beq");
    std::string phi_ic;
    int i = 0;
    while (getline(initial_condition, phi_ic)) {
        if (i >= n) {
            initial_condition.close();
            throw new std::runtime_error("Too many values in IC");
        }
        phi[i] = std::stod(phi_ic); 
        i++;
    }
    initial_condition.close();
    if (i != n){
        throw new std::runtime_error("Too few values in IC");
    }
}

int main() {
    print_header();


    // unpack some config data
    Config config = Config("config/config.json");
    Domain domain(config.domain_json());
    Equation * equation = create_equation(config.equation_json());
    Solver * solver = make_solver(config.solver_json(), domain);

    solver->set_initial_condition();
    solver->solve(*equation, domain);  

    delete solver;
    delete equation;
    return 0;
}
