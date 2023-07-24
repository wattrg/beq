#include <fstream>
#include "runge_kutta.h"

RungeKutta::RungeKutta(json json_config, Domain &domain, Equation *equation) 
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
    int n_comp = equation->number_components();
    unsigned n = domain.number_cells();
    unsigned n_ghost = domain.number_ghost();
    this->_phi_cpu = new Field<double>(n, n_ghost, n_comp, false); // allocated on CPU
    this->_residual = new Field<double>(n, n_ghost, n_comp, true); // allocated on GPU
    for (int stage = 0; stage < _n_stages; stage++) {
        Field<double> * field = new Field<double>(n, n_ghost, n_comp, true);
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
void apply_residual(double *phi, double *phi_new, double *residual, double dt, int nc, int nv){
    int index = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int stride = blockDim.x * gridDim.x;

    for (int ci = index; ci < nc + 1; ci += stride) {
        for (int vi = 0; vi < nv; vi++){
            int v_index = vi*(nc+2);
            phi_new[v_index+ci] =  phi[v_index+ci] + residual[v_index+ci] * dt;
        }
    }
}

StepResult RungeKutta::_take_step(Equation &equation, Domain &domain) {
    unsigned n_blocks = domain.number_blocks();
    unsigned block_size = domain.block_size();


    double dt = _cfl * equation.allowable_dt(*_phi_buffers[0], domain);
    for (int stage = 0; stage < _n_stages; stage++){
        fill_boundaries(*_phi_buffers[stage], domain, equation);

        equation.eval_residual(*_phi_buffers[stage], *_residual, domain);

        apply_residual<<<n_blocks,block_size>>>(
            _phi_buffers[stage]->data(), _phi_buffers[stage]->data(), _residual->data(), 
            dt, domain.number_cells(), equation.number_components()
        );
        auto code = cudaGetLastError();
        if (code != cudaSuccess) {
            std::cerr << "Cuda error in RungeKutta step: " << cudaGetErrorString(code) << std::endl;
            return StepResult::Failure;
        }

    }

    bool phi_valid = equation.check_phi(*_phi_buffers[0], domain);
    if (!phi_valid) {
        std::cerr << "Invalid state" << std::endl;
        return StepResult::Failure;
    }

    _t += dt;
    _time_since_last_plot += dt;

    return StepResult::Success;
}

bool RungeKutta::_stop() {
    return (_t >= _max_time || _step_number() >= _max_steps);
}

std::string RungeKutta::_stop_reason() {
    std::string reason;
    if (_t >= _max_time) {
        reason = std::string("Reached max_time");
    }
    else if (_step_number() >= _max_steps) {
        reason = std::string("Reached max_step");
    }
    else {
        reason = std::string("Unknown reason for stopping");
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
    unsigned start = _phi_cpu->number_components();
    unsigned stop = (_phi_cpu->length()+1) * _phi_cpu->number_components();
    for (unsigned i = start; i < stop; i++){
        file << std::setprecision(16) << (*_phi_cpu).flat_index(i);
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
    unsigned nc = _phi_cpu->length();
    unsigned nv = _phi_cpu->number_components();
    double *phi = _phi_cpu->data();

    // zero everything to begin.
    for (unsigned j = 0; j < _phi_cpu->size(); j++){
        phi[j] = 0.0;
    }

    unsigned count = 0;
    for (unsigned vi = 0; vi < _phi_cpu->number_components(); vi++){
        for (unsigned ci = 0; ci < _phi_cpu->length(); ci++) {
            if (getline(initial_condition, phi_ic)) {
                int index = vi*(nc+2) + ci + 1;
                phi[index] = std::stod(phi_ic);
                count++;
            }        
            else {
                std::cerr << "Too few values in IC" << std::endl;
                throw new std::runtime_error("Too few values in IC");
            }
        }
    }

    // while (getline(initial_condition, phi_ic)) {
    //     if (i > (_phi_cpu->length()+1)*_phi_cpu->number_components()) {
    //         initial_condition.close();
    //         std::cerr << "Too many values in IC" << std::endl;
    //         throw new std::runtime_error("Too many values in IC");
    //     }
    //     phi[i] = std::stod(phi_ic); 
    //     i++;
    // }
    initial_condition.close();
    if (count != nc*nv){
        std::cerr << "Incorrect number of values in IC" << std::endl;
        throw new std::runtime_error("Incorrect number of values in IC");
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
    for (unsigned j = 0; j < _phi_cpu->size(); j++) {
        std::cout << "phi[" << j << "] = " << phi[j] << "\n";
    }
    std::cout << std::endl;
}
