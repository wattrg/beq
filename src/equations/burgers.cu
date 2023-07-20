#include "equation.h"
#include <pthread.h>

Burgers::Burgers() {
    auto code = cudaMalloc(&_min_dt_gpu, sizeof(int));

    if (code != cudaSuccess){
        std::cerr << "Failed to allocate min_dt on GPU: " << cudaGetErrorString(code) << std::endl;
        throw std::runtime_error("Burgers: Could not allocate min_dt on GPU");
    }
}

Burgers::~Burgers() {
    auto code = cudaFree(_min_dt_gpu);

    if (code != cudaSuccess) {
        std::cerr << "Failed to free 'min_dt_gpu' in Burgers equation: "
                  << cudaGetErrorString(code)
                  << std::endl;
    }
}

__global__
void eval_burgers_residual(double *phi, double *residual, double dx, unsigned n) {
    int index = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int stride = blockDim.x * gridDim.x;

    for (unsigned i = index; i < n + 2; i += stride) {
        double flux_minus, flux_plus;
        flux_minus = phi[i-1] * phi[i-1];
        flux_plus = phi[i] * phi[i];
        residual[i] = - 0.5 * (flux_plus - flux_minus) / dx;
    }
}

void 
Burgers::eval_residual(Field<double> &phi, 
                       Field<double> &residual, 
                       Domain &domain) 
{
    unsigned n_blocks = domain.number_blocks();
    unsigned block_size = domain.block_size();

    eval_burgers_residual<<<n_blocks,block_size>>>(
        phi.data(), residual.data(), domain.dx(), phi.size()
    );

    auto code = cudaGetLastError();
    if (code != cudaSuccess) {
        std::cerr << "Cuda error in burgers residual eval: " 
                  << cudaGetErrorString(code) 
                  << std::endl;
        throw new std::runtime_error("Cuda error in burgers residual eval");
    }
}

__global__ 
void burgers_allowable_dt(double *phi, int n, int *min_dt, double dx) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    
    *min_dt = dx / phi[0] * 1000;
    for (int i = index; i < n; i +=  stride) {
        int dt = dx / phi[i] * 1000;
        atomicMin(min_dt, dt);
    }
}

double 
Burgers::allowable_dt(Field<double> &phi, Domain &domain) {
    unsigned n_blocks = domain.number_blocks();
    unsigned block_size = domain.block_size();

    burgers_allowable_dt<<<n_blocks,block_size>>>(
        phi.data(), phi.size(), _min_dt_gpu, domain.dx()
    );

    auto code = cudaGetLastError();
    if (code != cudaSuccess){
        std::cerr << "Burgers: Cuda error in 'allowable_dt': "
                  << cudaGetErrorString(code)
                  << std::endl;
        throw std::runtime_error("Cuda error in burgers allowable_dt");
    }

    code = cudaMemcpy(&_min_dt_cpu, _min_dt_gpu, sizeof(int), cudaMemcpyDeviceToHost);
    if (code != cudaSuccess){
        std::cerr << "Burgers: Cuda memcpy error in 'allowable_dt': "
                  << cudaGetErrorString(code)
                  << std::endl;
        throw std::runtime_error("Cuda error in burgers allowable_dt");
    }

    return (double)_min_dt_cpu / 1000.0;
}
