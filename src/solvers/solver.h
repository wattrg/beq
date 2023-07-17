#ifndef __SOLVER_H_
#define __SOLVER_H_

#include <iostream>
#include <string>

#include "../equations/equation.h"
#include "../domain.h"
#include "../boundary.h"

enum class StepResult {Success, Failure};

class Solver {
public:
    virtual ~Solver() {};
    void solve(Equation &equation, Domain &domain) {
        while (true) {
            StepResult step_result = _take_step(equation, domain);

            if (step_result == StepResult::Failure) {
                // for now, we'll just give up
                std::cerr << "Failed to take a step at step " << _step_num << std::endl;
                throw std::runtime_error("Failed to take a step");
            }

            _step_num++;

            if (_print_this_step()) {
                std::cout << "  " << _print_progress() << std::endl;
            }

            bool stop = _stop();
            if (_plot_this_step() || stop) {
                _write_solution();
            }

            if (stop) {
                std::cout << "Stopping: " << _stop_reason() << std::endl;
                break;
            }

        }
    }


    virtual void set_initial_condition() = 0;


protected:
    int _step_number() const { return _step_num; }
    virtual StepResult _take_step(Equation &equation, Domain &domain) = 0;
    virtual bool _stop() = 0;
    virtual std::string _stop_reason() = 0;
    virtual bool _print_this_step() = 0;
    virtual std::string _print_progress() = 0;
    virtual bool _plot_this_step() = 0;
    virtual void _write_solution() = 0;

private:
    int _step_num = 0;
};

Solver * make_solver(json solver_json, Domain &domain, Equation *equation);

#endif
