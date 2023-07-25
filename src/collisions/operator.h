#ifndef OPERATOR_H_
#define OPERATOR_H_

#include "../field.h"
#include "frequency.h"
#include "../domain.h"
#include "../equations/equation.h"

class CollisionOperator {
public:
    virtual ~CollisionOperator() {};

    virtual void 
    collide(Field<double> &phi, Domain &domain, Equation &equation) = 0;
};

class BGK : public CollisionOperator {
public:
    ~BGK();
    BGK(CollisionFrequency *frequency);

    void collide(Field<double> &phi, Domain &domain, Equation &eqation);

private:
    CollisionFrequency *_frequency_gpu;
    CollisionFrequency *_frequency_cpu;
};

#endif
