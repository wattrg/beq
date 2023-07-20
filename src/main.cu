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

#ifdef NDEBUG
#define DEBUG off
#else
#define DEBUG on
#endif

const char *LONG_HELP = 
    "USAGE:\n"
    "    beq command [options]\n"
    "\n"
    "Available commands:\n"
    "    - help: print this help message\n"
    "    - check: check the environment is set up correctly\n"
    "    - prep <CASE>: prepare a case\n"
    "    - run: run a case\n"
    "    - post: post process a case\n";

void print_header() {
    std::cout << "beq: Boltzmann equation solver\n";
    std::cout << "Git branch: " << STRINGIFY(GIT_BRANCH) << std::endl;
    std::cout << "Git commit: " << STRINGIFY(GIT_HASH) << std::endl;
    std::cout << "Build date: " << STRINGIFY(COMPILE_TIME) << std::endl;
    std::cout << "Debugging: " << STRINGIFY(DEBUG) << std::endl;
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

int check() {
    const char* beq = std::getenv("BEQ");
    if (!beq) {
        std::cerr << "Error: BEQ environment variable is not set" << std::endl;
        return 1;
    }

    return 0;
}

int prep(std::string case_name) {
    const char* beq = std::getenv("BEQ");
    if (!beq) {
        std::cerr << "Make sure BEQ environment variable is set" << std::endl;
        return 1;
    }
    std::string prep_name = std::string(beq) + "/lib/prep.py";
    std::string prep_command = "python " + prep_name + " " + case_name;
    int flag = std::system(prep_command.c_str());

    return flag;
}

int run() {
    // unpack some config data
    Config config = Config("config/config.json");

    Domain domain(config.domain_json());
    Equation * equation = create_equation(config.equation_json());
    Solver * solver = make_solver(config.solver_json(), domain, equation);

    solver->set_initial_condition();
    int result = solver->solve(*equation, domain);  

    delete solver;
    delete equation;

    return result;
}

int post(int argc, char* argv[]){
    const char* beq = std::getenv("BEQ");
    if (!beq) {
        std::cerr << "Make sure BEQ environment variable is set" << std::endl;
        return 1;
    }
    std::string post_command = "python " + std::string(beq) + "/lib/post.py";
    for (int i = 2; i < argc; i++){
        post_command.append(" ");
        post_command.append(argv[i]);
    }
    int flag = std::system(post_command.c_str());

    return flag;
}

enum class Command {Help, Check, Prep, Run, Post, Clean};

Command string_to_command(std::string command_string){
    if (command_string == "help") {
        return Command::Help;
    }
    if (command_string == "check") {
        return Command::Check;
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
    if (command_string == "clean") {
        return Command::Clean;
    }
    throw std::runtime_error("Unknown command");
}

int main(int argc, char* argv[]) {
    print_header();

    if (argc < 2){
        std::cout << std::endl;
        std::cout << LONG_HELP << std::endl;
        return 0;
    }

    Command command = string_to_command(std::string(argv[1]));
    int flag = 0;
    switch (command) {
        case Command::Help:
            std::cout << std::endl;
            std::cout << LONG_HELP << std::endl;
            break;
        case Command::Check:
            flag = check();
            if (flag == 0){
                std::cout << "Environment set up correctly" << std::endl;
            }
            break;
        case Command::Prep:
            std::cout << "Action: prepping " << std::string(argv[2]) << std::endl;
            flag = prep(std::string(argv[2]));
            break;
        case Command::Run:
            std::cout << "Action: running simulation" << std::endl;
            flag = run();
            break;
        case Command::Post:
            std::cout << "Action: post processing" << std::endl;
            flag = post(argc, argv);
            break;
        case Command::Clean:
            std::filesystem::remove_all("config");
            std::filesystem::remove_all("plot");
            std::filesystem::remove_all("solution");
    }
    return flag;
}
