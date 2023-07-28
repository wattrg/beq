#include "operator.h"
#include "../constants.h"
#include <json.hpp>
#include <stdexcept>

BGK::BGK(CollisionFrequency *frequency) {
    _frequency_cpu = frequency;
    auto code = cudaMalloc(&_frequency_gpu, sizeof(*frequency));
    if (code != cudaSuccess) {
        std::cerr << "cudaMalloc failure: " << cudaGetErrorString(code) << std::endl;
        throw new std::runtime_error("");
    }

    // assume that the CollisionFrequency model doesn't contain pointers else where
    code = cudaMemcpy(_frequency_gpu, _frequency_cpu, sizeof(*frequency), cudaMemcpyHostToDevice);
    if (code != cudaSuccess) {
        std::cerr << "cudaMemcpy failure: " << cudaGetErrorString(code) << std::endl;
        throw new std::runtime_error("");
    }
}

BGK::~BGK() {
    delete _frequency_cpu;
    cudaFree(_frequency_gpu);
}

__global__ void bgk_collide(double *phi, double *residual, int nc, int nv, double min_v, 
                            double dv, double mass, double volume, CollisionFrequency *frequency, double r) {
    int index = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int stride = blockDim.x * gridDim.x;


    double R = kB / mass;

    for (int ci = index; ci < nc + 1; ci += stride) {
        double n_particles = 0.0;
        double v_avg = 0.0;
        double thermal_energy = 0.0;

        // step 1: Compute moments of the distribution for this cell
        // compute number of particles and velocity
        for (int vi = 0; vi < nv; vi++) {
            double v = min_v + (vi + 0.5) * dv;
            int v_index = vi*(nc+2);

            n_particles += phi[v_index + ci] * dv;
            v_avg += v * phi[v_index + ci] * dv; 
        } 
        v_avg /= n_particles;

        // thermal energy
        for (int vi = 0; vi < nv; vi++) {
            double v = min_v + (vi + 0.5) * dv;
            int v_index = vi*(nc+2);

            thermal_energy += (v - v_avg) * (v - v_avg) * phi[v_index + ci] * dv;
        }
        thermal_energy = mass * thermal_energy / (2 * volume);

        // step 2: macroscopic properties for this cell
        double density = mass * n_particles / volume;
        double num_density = n_particles / volume;
        double temp = thermal_energy / (density * 0.5 * R);
        // printf("v_avg = %.16e, density = %.16e, temp = %.16e\n", v_avg, density, temp);

        // step 3: collision frequency for this cell
        double mu = 16.0 * num_density * r * r * sqrt(PI * R * temp);
        // double mu = frequency->collision_frequency(temp, n_particles);

        // step 4: apply collisions for this cell
        for (int vi = 0; vi < nv; vi++) {
            double v = min_v + (vi + 0.5) * dv;
            int v_index = vi*(nc+2);

            // compute equilibrium distribution value
            double exponent = -mass / (2*kB*temp) * (v-v_avg) * (v-v_avg); 
            double norm = sqrt(mass / (2 * PI * kB * temp));
            double f0 = num_density * norm * exp(exponent);

            // compute the effect of collisions
            residual[v_index+ci] += mu * num_density * (f0 - phi[v_index+ci]);
            printf("mu = %.16e, f0 = %.16e, f = %.16e residual = %.16e\n", mu, f0, phi[v_index+ci], f0 - phi[v_index+ci]);
        }
    }
}

void BGK::collide(Field<double> &phi, Field<double> &residual, 
                  Domain &domain, Equation &equation, double min_v, 
                  double dv, double mass, double r) 
{
    unsigned n_blocks = domain.number_blocks();
    unsigned block_size = domain.block_size();

    bgk_collide<<<n_blocks, block_size>>>(
        phi.data(), residual.data(), domain.number_cells(), 
        equation.number_components(), min_v, dv, mass, 
        domain.dx(), _frequency_gpu, r
    );

    auto code = cudaGetLastError();
    if (code != cudaSuccess) {
        std::cerr << "Cuda error in BGK collision term: " << cudaGetErrorString(code) << std::endl;
        throw new std::runtime_error("Cuda error");
    }
}
