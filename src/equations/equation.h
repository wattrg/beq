#ifndef EQUATION_H_
#define EQUATION_H_

#include <json.hpp>
#include "../field.h"   
#include "../domain.h"
#include "../collisions/operator.h"

using json = nlohmann::json;

class Equation {
public:
    virtual ~Equation() {}
    virtual void eval_residual(Field<double> &phi, Field<double> &residual, Domain &domain) = 0;
    virtual double allowable_dt(Field<double> &phi, Domain &domain) = 0;
    virtual int number_components() = 0;        
    virtual bool check_phi(Field<double> &phi, Domain &domain) { (void)phi; (void)domain; return true; };
};

class Advection : public Equation {
public:
    Advection(double velocity);
    Advection(json json_data);

    void eval_residual(Field<double> &phi, Field<double> &residual, Domain &domain); 

    double allowable_dt(Field<double> &phi, Domain &domain) { 
        // trick the compiler into thinking we used phi, so we don't
        // get unused parameter warnings when compiling
        (void)phi;

        // compute the allowable time step
        return domain.dx() / _velocity; 
    }

    int number_components() {return 1;}

private:
    double _velocity;
};

class Burgers : public Equation {
public:
    Burgers();
    ~Burgers();

    void eval_residual(Field<double> &phi, Field<double> &residual, Domain &domain);

    double allowable_dt(Field<double> &phi, Domain &domain);

    int number_components() {return 1;}

private:
    int *_min_dt_gpu;
    int _min_dt_cpu;
};

class CollisionOperator;

class Boltzmann : public Equation {
public:
    Boltzmann(json json_data, json gas_model);
    ~Boltzmann(){
        cudaFree(_phi_valid_gpu);
    }

    void eval_residual(Field<double> &phi, Field<double> &residual, 
                       Domain &domain);

    double allowable_dt(Field<double> &phi, Domain &domain);

    int number_components() {
        return _nv;
    }

    bool check_phi(Field<double> &phi, Domain &domain);

private:
    CollisionOperator *_collisions;
    double _mass;
    double _r;
    double _dv;
    double _min_v;
    double _max_v;
    int _nv;
    bool *_phi_valid_gpu;
};

Equation * create_equation(json json_data, json gas_model);

#endif
