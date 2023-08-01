#include <iostream>
#include "equation.h"

Diffusion::Diffusion(double diffusivity) : _diffusivity(diffusivity) {}

Diffusion::Diffusion(json json_data) {
    _diffusivity = json_data.at("diffusivity");
}
__global__
void eval_diffusion_residual(double *phi, double *residual, double k, double dx, int n) {
    int index = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int stride = blockDim.x * gridDim.x;

    for (int i = index; i < n+1; i += stride) {
        residual[i] = k * (phi[i+1] - 2*phi[i] + phi[i-1]) / (dx*dx);
    }
}

void Diffusion::eval_residual(Field<double> &phi, Field<double> &residual, Domain &domain) {
    unsigned n_blocks = domain.number_blocks();
    unsigned block_size = domain.block_size();

    eval_diffusion_residual<<<n_blocks,block_size>>>(
        phi.data(), residual.data(), _diffusivity, domain.dx(), phi.size()
    );

    auto code = cudaGetLastError();
    if (code != cudaSuccess){
        std::cerr << "Cuda error in diffusion residual eval: " << cudaGetErrorString(code) << std::endl;
        throw new std::runtime_error("Encountered cuda error");
    }
}
