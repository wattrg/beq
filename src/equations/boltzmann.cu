#include "equation.h"

Boltzmann::Boltzmann(json json_data) {
    _dv = json_data.at("dv");
    _n_velocity_increments = json_data.at("n_vel_increments");
    _max_v = _n_velocity_increments * _dv;
}

__global__
void eval_boltzmann_residual(double *phi, double *residual, 
                             double dx, double dv, int n, int nv)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    for (int xi = index; xi < n; xi += stride) {
        for (int vi = 0; vi < nv; vi++) {
            double phi_minus, phi_plus;

            // right moving particles
            if (xi == 0) {
                phi_minus = phi[vi*n + n-1];
                phi_plus = phi[vi*n + xi];
            }
            else {
                phi_minus = phi[vi*n + xi - 1];
                phi_plus = phi[vi*n + xi];
            }
            residual[vi*n + xi] = (vi + 0.5) * dv * (phi_plus - phi_minus);

            // left moving particles
            if (xi == n-1) {
                phi_plus = phi[(vi+nv)*n + 0];
                phi_minus = phi[(vi+nv)*n + xi];
            }
            else {
                phi_minus = phi[(vi+nv)*n + xi];
                phi_plus = phi[(vi+nv)*n + xi + 1];
            }
            residual[(vi+nv)*n + xi] = -(vi + 0.5) * dv * (phi_plus - phi_minus);
        } 
    }
}

void Boltzmann::eval_residual(Field<double> &phi, Field<double> &residual, 
                              Domain &domain) 
{
    unsigned n_blocks = domain.number_blocks();
    unsigned block_size = domain.block_size();

    eval_boltzmann_residual<<<n_blocks, block_size>>>(
        phi.data(), residual.data(), domain.dx(), _dv, 
        domain.number_cells(), _n_velocity_increments
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
