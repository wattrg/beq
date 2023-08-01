#include "equation.h"

Equation * create_equation(json json_data, json gas_model){
    Equation * eq;
    std::string type = json_data.at("type");
    if (type == "advection"){
        eq = new Advection(json_data);
    }
    else if (type == "diffusion") {
        eq = new  Diffusion(json_data);
    }
    else if (type == "burgers") {
        eq = new Burgers();
    }
    else if (type == "boltzmann") {
        eq = new Boltzmann(json_data, gas_model);
    }
    else {
        throw new std::runtime_error("Unknwon equation type");
    }
    return eq;
}
