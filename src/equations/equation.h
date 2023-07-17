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
    virtual double allowable_dt(Field<double> &phi, Domain &domain) = 0;
    virtual int number_components() = 0;        
    virtual bool check_phi(Field<double> &phi, Domain &domain) = 0;
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

    bool check_phi(Field<double> &phi, Domain &domain) { return true; }

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

    bool check_phi(Field<double> &phi, Domain &domain) { return true; }

private:
    int *_min_dt_gpu;
    int _min_dt_cpu;
};

class Boltzmann : public Equation {
public:
    Boltzmann(json json_data);
    ~Boltzmann(){}

    void eval_residual(Field<double> &phi, Field<double> &residual, Domain &domain);

    double allowable_dt(Field<double> &phi, Domain &domain);

    int number_components() {
        return _nv;
    }

    bool check_phi(Field<double> &phi, Domain &domain);

private:
    double _dv;
    double _min_v;
    double _max_v;
    int _nv;
    bool *_phi_valid_gpu;
};

Equation * create_equation(json json_data);

#endif
