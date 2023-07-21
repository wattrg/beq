#include <iostream>
#include "equation.h"

Advection::Advection(double velocity) : _velocity(velocity) {}

Advection::Advection(json json_data) {
    _velocity = json_data.at("velocity");
}
__global__
void eval_advection_residual(double *phi, double *residual, double u, double dx, int n) {
    int index = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int stride = blockDim.x * gridDim.x;

    for (int i = index; i < n+2; i += stride) {
        double phi_minus, phi_plus;
        if (u > 0) {
            phi_plus = phi[i];
            phi_minus = phi[i-1];
        }
        else {
            phi_plus = phi[i+1];
            phi_minus = phi[i];
        }
        residual[i] = -u * (phi_plus - phi_minus) / dx;
    }
}

void Advection::eval_residual(Field<double> &phi, Field<double> &residual, Domain &domain) {
    unsigned n_blocks = domain.number_blocks();
    unsigned block_size = domain.block_size();

    eval_advection_residual<<<n_blocks,block_size>>>(
        phi.data(), residual.data(), _velocity, domain.dx(), phi.size()
    );

    auto code = cudaGetLastError();
    if (code != cudaSuccess){
        std::cerr << "Cuda error in advection residual eval: " << cudaGetErrorString(code) << std::endl;
        throw new std::runtime_error("Encountered cuda error");
    }
}
