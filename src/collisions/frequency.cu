#include "frequency.h"
#include <stdexcept>
#include <iostream>


CollisionFrequency * make_collision_frequency(json frequency, json gas_model) {
    std::string type = frequency.at("type");
    if (type == "hard_sphere") {
        double radius = gas_model.at("radius");
        double mass = gas_model.at("mass");
        return new HardSphere(radius, mass);
    }
    std::cerr << "Unknown collision frequency: " << type << std::endl;
    throw std::runtime_error("Unknown collision frequency");
}
