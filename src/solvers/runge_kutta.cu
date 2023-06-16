#include <fstream>
#include "runge_kutta.h"

RungeKutta::RungeKutta(json json_config, Domain &domain) 
    : _t(0.0), _time_since_last_plot(0.0)
{
    // read configuration
    this->_cfl = json_config.at("cfl"); 
    this->_n_stages = json_config.at("number_stages");
    this->_max_time = json_config.at("max_time");
    this->_max_steps = json_config.at("max_step");
    this->_print_frequency = json_config.at("print_frequency");
    this->_plot_every_n_steps = json_config.at("plot_every_n_steps");
    this->_plot_frequency = json_config.at("plot_frequency");

    // allocate memory
    unsigned n = domain.number_cells();
    this->_phi_cpu = new Field<double>(n, false); // allocated on CPU
    this->_residual = new Field<double>(n, true); // allocated on GPU
    for (int stage = 0; stage < _n_stages; stage++) {
        Field<double> * field = new Field<double>(n, true);
        this->_phi_buffers.push_back(field);
    }
}


RungeKutta::~RungeKutta() {
    delete this->_residual;
    delete this->_phi_cpu;
    for (int stage = 0; stage < _n_stages; stage++){
        delete this->_phi_buffers[stage];
    }
}

__global__
void apply_residual(double *phi, double *phi_new, double *residual, double dt, int n){
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    for (int i = index; i < n; i += stride) {
        phi_new[i] =  phi[i] + residual[i] * dt;
    }
}

void RungeKutta::_take_step(Equation &equation, Domain &domain) {
    unsigned n_blocks = domain.number_blocks();
    unsigned block_size = domain.block_size();

    double dt = _cfl * equation.allowable_dt(*_phi_buffers[0], domain);

    for (int stage = 0; stage < _n_stages; stage++){
        equation.eval_residual(*_phi_buffers[stage], *_residual, domain);

        apply_residual<<<n_blocks,block_size>>>(
            _phi_buffers[stage]->data(), _phi_buffers[stage]->data(), _residual->data(), 
            dt, domain.number_cells()
        );
        auto code = cudaGetLastError();
        if (code != cudaSuccess) {
            std::cerr << "Cuda error in RungeKutta step: " << cudaGetErrorString(code) << std::endl;
            throw new std::runtime_error("Encountered cuda error");
        }

    }


    _t += dt;
    _time_since_last_plot += dt;
}

bool RungeKutta::_stop() {
    return (_t >= _max_time || _step_number() >= _max_steps);
}

std::string RungeKutta::_stop_reason() {
    std::string reason;
    if (_t >= _max_time) {
        reason = std::string("Reached max_time");
    }
    if (_step_number() >= _max_steps) {
        reason = std::string("Reached max_step");
    }
    return reason;
}

bool RungeKutta::_print_this_step() {
    return (_step_number() % _print_frequency == 0);
}

void RungeKutta::_write_solution() {
    // copy solution to CPU
    cudaMemcpy(_phi_cpu->data(), _phi_buffers[0]->data(), _phi_cpu->memory(), cudaMemcpyDeviceToHost);

    // write contents of cpu buffer to file
    std::string file_name = "solution/phi_" + std::to_string(_n_solutions) + ".beq";
    std::ofstream file(file_name);
    for (unsigned i = 0; i < _phi_cpu->size(); i++){
        file << (*_phi_cpu)(i);
        file << "\n";
    }
    file.close();
    _time_since_last_plot = _t - _n_solutions * _plot_frequency;
    _n_solutions++;
}

std::string RungeKutta::_print_progress() {
    std::string progress ("Step: ");
    progress.append(std::to_string(_step_number()));
    progress.append(", t = ");
    progress.append(std::to_string(_t));
    return progress;
}

bool RungeKutta::_plot_this_step() {
    bool plot_step = _plot_every_n_steps > 0 && _step_number() % _plot_every_n_steps == 0;
    bool plot_time = _plot_frequency > 0 && _time_since_last_plot >= _plot_frequency;
    return (plot_step || plot_time);
}

void RungeKutta::set_initial_condition() {
    // read initial condition from file
    std::ifstream initial_condition("solution/phi_0.beq");
    std::string phi_ic;
    unsigned i = 0;
    double *phi = _phi_cpu->data();
    while (getline(initial_condition, phi_ic)) {
        if (i >= _phi_cpu->size()) {
            initial_condition.close();
            throw new std::runtime_error("Too many values in IC");
        }
        phi[i] = std::stod(phi_ic); 
        i++;
    }
    initial_condition.close();
    if (i != _phi_cpu->size()){
        throw new std::runtime_error("Too few values in IC");
    }

    // copy initial condition to GPU buffer
    auto code = cudaMemcpy(
        _phi_buffers[0]->data(), _phi_cpu->data(), _phi_cpu->memory(), 
        cudaMemcpyHostToDevice
    );

    if (code != cudaSuccess) {
        std::cerr << "Cuda error setting initial condition: " << cudaGetErrorString(code) << std::endl;
        throw new std::runtime_error("Encountered cuda error");
    }
}
