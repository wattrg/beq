#include "equation.h"

Boltzmann::Boltzmann(json json_data) {
    _min_v = json_data.at("min_v");
    _max_v = json_data.at("max_v");
    _nv = json_data.at("n_vel_increments");
    _dv = (_max_v - _min_v) / _nv;
}

__global__
void eval_boltzmann_residual(double *phi, double *residual, 
                             double dx, double dv, int n, int nv, double min_v)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    for (int xi = index; xi < n; xi += stride) {
        for (int vi = 0; vi < nv; vi++) {
            double phi_minus, phi_plus;
            double v = min_v + (vi + 0.5) * dv;

            // left boundary
            if (xi == 0) {
                if (v > 0) {
                    // right moving particles
                    phi_minus = phi[vi*n + n-1];
                    phi_plus = phi[vi*n + xi];
                }
                else{
                    // left moving particles
                    phi_plus = phi[vi*n + xi + 1];
                    phi_minus = phi[vi*n+ xi];
                }
            }

            // right boundary
            if (xi == n-1) {
                if (v < 0) {
                    // left moving particle
                    phi_plus = phi[vi*n + 0];
                    phi_minus = phi[vi*n + xi];
                }
                else {
                    // right moving particle
                    phi_plus = phi[vi*n + xi];
                    phi_minus = phi[vi*n +xi - 1];
                }
            }

            // interior cell
            if (v < 0) {
                // left moving particles
                phi_plus = phi[vi*n + xi + 1];
                phi_minus = phi[vi*n + xi];
            }
            else {
                // right moving particles
                phi_plus = phi[vi*n + xi];
                phi_minus = phi[vi*n + xi - 1];
            }

            residual[vi*n + xi] = v * (phi_plus - phi_minus);
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
}

double Boltzmann::allowable_dt(Field<double> &phi, Domain &domain){
    (void) phi;
    return domain.dx() / _max_v;
}
