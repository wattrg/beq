#include <iostream>
#include <fstream>
#include <math.h>
#include <filesystem>

#include "config.h"

#define _STRINGIFY(x) #x
#define STRINGIFY(x) _STRINGIFY(x)

const int BLOCK_SIZE = 256;
const double PI = 3.14159265;
int number_solutions = 0;

void print_header() {
    std::cout << "beq: Boltzmann equation solver\n";
    std::cout << "Git branch: " << STRINGIFY(GIT_BRANCH) << "\n";
    std::cout << "Git commit: " << STRINGIFY(GIT_HASH) << "\n";
    std::cout << "Build date: " << STRINGIFY(COMPILE_TIME) << "\n";
}

__global__
void eval_rhs(double *phi, double *residual, double u, double dx, int n) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    for (int i = index; i < n; i += stride) {
        double phi_minus, phi_plus;
        if (i == 0) {
            phi_minus = phi[n-1];
            phi_plus = phi[i];
        }
        else {
            phi_minus = phi[i-1];
            phi_plus = phi[i];
        }
        residual[i] = -u * (phi_plus - phi_minus) / dx;
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

void take_step(double *phi, double *phi_new, double *residual, double u, double dx, double dt, int n, int number_blocks) {
    eval_rhs<<<number_blocks, BLOCK_SIZE>>>(phi, residual, u, dx, n);
    apply_residual<<<number_blocks, BLOCK_SIZE>>>(phi, phi_new, residual, dt, n);

    // swap phi and phi_new
    double *phi_tmp = phi;
    phi = phi_new;
    phi_new = phi_tmp;
}

void fill_initial_condition(double *phi, double dx, int n) {
    double sigma = 0.4;
    double x0 = 2.0;
    // initialise the data
    for (int i = 0; i < n; i++){
        double x = (i + 0.5) * dx;
        double exponent = 0.5 * (x - x0) / sigma;
        phi[i] = 1.0 / (sigma * sqrt(2*PI)) * exp(-exponent*exponent);
    }
}

void write_solution(double *phi_gpu, double *phi_cpu, int n){
    // copy data to CPU
    cudaMemcpy(phi_cpu, phi_gpu, n*sizeof(double), cudaMemcpyDeviceToHost);

    // write contents of cpu buffer to file
    std::string file_name = "solution/phi_" + std::to_string(number_solutions) + ".beq";
    std::ofstream file(file_name);
    for (int i = 0; i < n; i++){
        file << phi_cpu[i];
        file << "\n";
    }
    file.close();
    number_solutions++;
}

int main() {
    print_header();


    // unpack some config data
    Config config = Config("config/config.json");
    int N = config.number_cells();
    double L = config.length();
    double u = config.velocity();
    const int number_blocks = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    double dx = L / N;
    double dt = config.cfl() * dx / u;
    int print_frequency = config.print_frequency();
    int plot_frequency = config.plot_frequency();
    double max_time = config.max_time();
    int max_step = config.max_step();

    // allocate memory on host
    double *phi = new double[N];
    fill_initial_condition(phi, dx, N);

    // allocate memory on GPU
    double *phi_gpu, *phi_new_gpu, *residual_gpu;
    cudaMalloc(&phi_gpu, N*sizeof(double));
    cudaMalloc(&phi_new_gpu, N*sizeof(double));
    cudaMalloc(&residual_gpu, N*sizeof(double));


    // copy initial condition to the GPU
    cudaMemcpy(phi_gpu, phi, N*sizeof(double), cudaMemcpyHostToDevice);

    // write initial condition
    std::filesystem::create_directory("solution");
    write_solution(phi_gpu, phi, N);


    double t = 0;
    for (int step = 0; step < max_step; step++){
        take_step(phi_gpu, phi_gpu, residual_gpu, u, dx, dt, N, number_blocks);
        t += dt;

        if (step % plot_frequency == 0){
            write_solution(phi_gpu, phi, N);
        }

        if (step % print_frequency == 0){
            std::cout << "step: " << step << "\n";
        }

        if (t > max_time){
            std::cout << "Finished: reached max_time" << std::endl;
            break;
        }
    }
    write_solution(phi_gpu, phi, N);


    // free host memory
    delete [] phi;

    // free device memory
    cudaFree(phi_gpu);
    cudaFree(phi_new_gpu);
    cudaFree(residual_gpu);
    
    return 0;
}
