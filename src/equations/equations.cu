#include "equation.h"

Equation * create_equation(json json_data){
    Equation * eq;
    std::string type = json_data.at("type");
    if (type == "advection"){
        eq = new Advection(json_data);
    }
    else if (type == "burgers") {
        eq = new Burgers();
    }
    else {
        throw new std::runtime_error("Unknwon equation type");
    }
    return eq;
}
