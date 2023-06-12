#ifndef __EQUATION_H_
#define __EQUATION_H_

#include <json.hpp>
#include "../field.h"   
#include "../domain.h"

using json = nlohmann::json;

class Equation {
public:
    virtual ~Equation() {}
    virtual void eval_residual(Field<double> &phi, Field<double> &residual, Domain &domain) = 0;
    virtual double allowable_dt(Domain &domain) = 0;
        
};

class Advection : public Equation {
public:
    Advection(double velocity);
    Advection(json json_data);

    void eval_residual(Field<double> &phi, Field<double> &residual, Domain &domain); 
    double allowable_dt(Domain &domain) { return domain.dx() / _velocity; }

private:
    double _velocity;
};

Equation * create_equation(json json_data);

#endif
