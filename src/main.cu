#include <cstdlib>
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

const char *LONG_HELP = 
    "USAGE: "
    "    beq command [options]\n"
    "\n"
    "Available commands:\n"
    "    - help: print this help message\n"
    "    - prep <CASE>: prepare a case\n"
    "    - run: run a case\n"
    "    - post: post process a case\n";

void print_header() {
    std::cout << "beq: Boltzmann equation solver\n";
    std::cout << "Git branch: " << STRINGIFY(GIT_BRANCH) << std::endl;
    std::cout << "Git commit: " << STRINGIFY(GIT_HASH) << std::endl;
    std::cout << "Build date: " << STRINGIFY(COMPILE_TIME) << std::endl;
    std::cout << std::endl;
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

int prep(std::string case_name) {
    std::cout << "Preparing " << case_name << std::endl;
    const char* beq = std::getenv("BEQ");
    if (!beq) {
        std::cerr << "Make sure BEQ environment variable is set" << std::endl;
        return 1;
    }
    std::string prep_name = std::string(beq) + "/bin/beq_prep";
    std::string prep_command = "python " + prep_name + " " + case_name;
    std::system(prep_command.c_str());

    return 0;
}

int run() {
    // unpack some config data
    Config config = Config("config/config.json");
    std::cout << "Running " << config.title() << std::endl;

    Domain domain(config.domain_json());
    Equation * equation = create_equation(config.equation_json());
    Solver * solver = make_solver(config.solver_json(), domain);

    solver->set_initial_condition();
    solver->solve(*equation, domain);  

    delete solver;
    delete equation;

    return 0;
}

int post(){
    std::cout << "Post Processing" << std::endl;
    const char* beq = std::getenv("BEQ");
    if (!beq) {
        std::cerr << "Make sure BEQ environment variable is set" << std::endl;
        return 1;
    }
    std::string post_command = "python " + std::string(beq) + "/bin/beq_post";
    std::system(post_command.c_str());

    return 0;
}

enum class Command {Help, Prep, Run, Post};

Command string_to_command(std::string command_string){
    if (command_string == "help") {
        return Command::Help;
    }
    if (command_string == "prep") {
        return Command::Prep;
    }
    if (command_string == "run") {
        return Command::Run;
    }
    if (command_string == "post") {
        return Command::Post;
    }
    throw std::runtime_error("Unknown command");
}

int main(int argc, char* argv[]) {
    print_header();

    if (argc < 2){
        std::cout << LONG_HELP << std::endl;
        return 0;
    }

    Command command = string_to_command(std::string(argv[1]));
    switch (command) {
        case Command::Help:
            std::cout << LONG_HELP << std::endl;
            break;
        case Command::Prep:
            return prep(std::string(argv[2]));
            break;
        case Command::Run:
            return run();
            break;
        case Command::Post:
            return post();
            break;
    }
}
