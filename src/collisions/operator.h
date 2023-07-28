#ifndef OPERATOR_H_
#define OPERATOR_H_

#include "../field.h"
#include "frequency.h"
#include "../domain.h"
#include "../equations/equation.h"

#include <json.hpp>

using json = nlohmann::json;

class CollisionOperator {
public:
    virtual ~CollisionOperator() {};

    virtual void 
    collide(Field<double> &phi, Field<double> &residual, Domain &domain, 
            Equation &equation, double min_v, double dv, double mass, double r) = 0;
};

class BGK : public CollisionOperator {
public:
    ~BGK();
    BGK(CollisionFrequency *frequency);

    void collide(Field<double> &phi, Field<double> &residual, Domain &domain, 
                 Equation &eqation, double min_v, double dv, double mass, double r);

private:
    CollisionFrequency *_frequency_gpu;
    CollisionFrequency *_frequency_cpu;
};

CollisionOperator * make_collision_operator(json coll_op, json gas_model);

#endif
