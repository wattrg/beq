#ifndef __RUNGE_KUTTA_H_
#define __RUNGE_KUTTA_H_

#include <json.hpp>
#include <vector>
#include "../config.h"
#include "../field.h"
#include "solver.h"
#include "../domain.h"
#include "../equations/equation.h"

using json = nlohmann::json;

class RungeKutta : public Solver {
public:
    RungeKutta (json solver_config, Domain &domain);
    ~RungeKutta();

    void set_initial_condition();


private:
    int _n_stages;
    double _t;
    double _max_time;
    int _max_steps;
    int _print_frequency;
    int _plot_every_n_steps;
    double _plot_frequency;
    double _time_since_last_plot;
    int _n_solutions = 1;
    double _cfl;
    Field<double> *_phi_cpu;
    Field<double> *_residual;
    std::vector<Field<double>*> _phi_buffers;


    void _take_step(Equation &equation, Domain &domain);
    bool _stop();
    std::string _stop_reason();
    bool _print_this_step();
    std::string _print_progress();
    bool _plot_this_step();
    void _write_solution();
};

#endif

