#include "equation.h"

Boltzmann::Boltzmann(json json_data, json gas_model) {
    _min_v = json_data.at("min_v");
    _max_v = json_data.at("max_v");
    _nv = json_data.at("n_vel_increments");
    _dv = (_max_v - _min_v) / _nv;
    _mass = gas_model.at("mass");

    auto code = cudaMalloc(&_phi_valid_gpu, sizeof(bool));

    if (code != cudaSuccess) {
        std::cerr << "Failed to allocate phi_valid on the GPU: " << cudaGetErrorString(code) << std::endl;
        throw std::runtime_error("Failed to allocate phi_valid on GPU.");
    }
    bool valid = true;
    code = cudaMemcpy(_phi_valid_gpu, &valid, sizeof(valid), cudaMemcpyHostToHost);    
    if (code != cudaSuccess) {
        std::cerr << "Failed to set phi_valid to true: " << cudaGetErrorString(code) << std::endl;
        throw std::runtime_error("Failed to set phi_vaild to true");
    }

    auto collision_operator_json = json_data.at("collision_operator");
    _collisions = make_collision_operator(collision_operator_json, gas_model);
}

__global__
void eval_boltzmann_residual(double *phi, double *residual, 
                             double dx, double dv, int nc, int nv, double min_v)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int stride = blockDim.x * gridDim.x;

    for (int ci = index; ci < nc + 1; ci += stride) {
        for (int vi = 0; vi < nv; vi++) {
            double phi_minus, phi_plus;
            double v = min_v + (vi + 0.5) * dv;

            int v_index = vi*(nc+2);

            if (v < 0.0) {
                // left moving particles
                phi_plus = phi[v_index + ci + 1];
                phi_minus = phi[v_index + ci];
            }
            else if (v > 0.0) {
                // right moving particles
                phi_plus = phi[v_index + ci];
                phi_minus = phi[v_index + ci - 1];
            }
            else {
                // stationary particles
                phi_plus = phi[v_index + ci];
                phi_minus = phi[v_index + ci];
            }

            residual[v_index + ci] = - v * (phi_plus - phi_minus) / dx;
        } 
    }
}

void Boltzmann::eval_residual(Field<double> &phi, Field<double> &residual, 
                              Domain &domain) 
{
    unsigned n_blocks = domain.number_blocks();
    unsigned block_size = domain.block_size();

    eval_boltzmann_residual<<<n_blocks, block_size>>>(
        phi.data(), residual.data(), domain.dx(), _dv, domain.number_cells(), _nv, _min_v
    );

    auto code = cudaGetLastError();
    if (code != cudaSuccess) {
        std::cerr << "Cuda error in Boltzmann residual eval: " << cudaGetErrorString(code) << std::endl;
        throw new std::runtime_error("Encountered cuda error");
    }

    _collisions->collide(phi, residual, domain, *this, _min_v, _dv, _mass);
}

double Boltzmann::allowable_dt(Field<double> &phi, Domain &domain){
    (void) phi;
    double max_v = fmax(fabs(_min_v), fabs(_max_v));
    return domain.dx() / max_v;
}

__global__
void check_phi_gpu(double *phi, bool *valid, int nc, int nv) {
    int index = blockIdx.x * blockDim.x + threadIdx.x+1;
    int stride = blockDim.x * gridDim.x;

    for (int ci = index; ci < nc+1; ci += stride) {
        for (int vi = 0; vi < nv; vi++){
            if (phi[vi*(nc+2) + ci] < -1.0) {
                printf("ci = %d, vi = %d, phi = %.16e\n", ci, vi, phi[vi*nc + ci]);
                *valid = false;
            }
        }
    }
}

bool Boltzmann::check_phi(Field<double> &phi, Domain &domain) {
    unsigned n_blocks = domain.number_blocks();
    unsigned block_size = domain.block_size();

    check_phi_gpu<<<n_blocks, block_size>>>(phi.data(), _phi_valid_gpu, domain.number_cells(), _nv);

    auto code = cudaGetLastError();
    if (code != cudaSuccess) {
        std::cerr << "Cuda error while checking phi: " << cudaGetErrorString(code) << std::endl;
        return false;
    }

    bool phi_valid;
    code = cudaMemcpy(&phi_valid, _phi_valid_gpu, sizeof(bool), cudaMemcpyDeviceToHost);
    if (code != cudaSuccess) {
        std::cerr << "Boltzmann: Cuda memcpy error for phi_valid: " << cudaGetErrorString(code) << std::endl;
        return false;
    }

    return phi_valid;
}
