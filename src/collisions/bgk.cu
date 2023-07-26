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
                            double dv, double mass, double volume, CollisionFrequency *frequency) {
    int index = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int stride = blockDim.x * gridDim.x;


    double R = kB / mass;

    for (int ci = index; ci < nc + 1; ci += stride) {
        double n_particles = 0.0;
        double v_avg = 0.0;
        double thermal_energy = 0.0;

        // compute number of particles and velocity for this cell
        for (int vi = 0; vi < nv; vi++) {
            double v = min_v + (vi + 0.5) * dv;
            int v_index = vi*(nc+2);

            n_particles += phi[v_index + ci];
            v_avg += v * phi[v_index + ci]; 
        } 
        v_avg /= n_particles;

        // compute thermal energy for this cell
        for (int vi = 0; vi < nv; vi++) {
            double v = min_v + (vi + 0.5) * dv;
            int v_index = vi*(nc+2);

            thermal_energy += (v - v_avg) * (v - v_avg) * phi[v_index + ci];
        }
        thermal_energy /= n_particles;

        // macroscopic properties for this cell
        double density = mass * n_particles / volume;
        double temp = 2 * thermal_energy / (density * R);

        // collision frequency for this cell
        printf("Before calling collision frequency");
        double mu = frequency->collision_frequency(temp, n_particles);
        printf("After calling collision frequency");

        // apply collisions for this cell
        for (int vi = 0; vi < nv; vi++) {
            double v = min_v + (vi + 0.5) * dv;
            int v_index = vi*(nc+2);

            // compute equilibrium distribution value
            double exponent = -mass / (2*kB*temp) * (v-v_avg) * (v-v_avg); 
            double norm = sqrt(mass / (2 * PI * kB * temp));
            double f0 = n_particles * norm * exp(exponent);

            residual[v_index+ci] += n_particles * mu * (f0 - phi[v_index+ci]);
        }
    }
}

void BGK::collide(Field<double> &phi, Field<double> &residual, 
                  Domain &domain, Equation &equation, double min_v, 
                  double dv, double mass) 
{
    unsigned n_blocks = domain.number_blocks();
    unsigned block_size = domain.block_size();

    bgk_collide<<<n_blocks, block_size>>>(
        phi.data(), residual.data(), domain.number_cells(), 
        equation.number_components(), min_v, dv, mass, 
        domain.dx(), _frequency_gpu
    );

    auto code = cudaGetLastError();
    if (code != cudaSuccess) {
        std::cerr << "Cuda error in BGK collision term: " << cudaGetErrorString(code) << std::endl;
        throw new std::runtime_error("Cuda error");
    }
}
