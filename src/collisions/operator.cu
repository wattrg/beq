#include "operator.h"

CollisionOperator * make_collision_operator(json coll_op, json gas_model) {
    std::string type = coll_op.at("type");
    if (type == "BGK") {
        json collision_frequency = coll_op.at("collision_frequency");
        CollisionFrequency *mu = make_collision_frequency(collision_frequency, gas_model);
        return new BGK(mu);
    }

    std::cerr << "Unknown collision operator: " << type << std::endl;
    throw new std::runtime_error("Unknown collision operator");
}
