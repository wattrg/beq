#include "operator.h"
#include "../constants.h"
#include <json.hpp>

BGK::BGK(CollisionFrequency *frequency) {
    _frequency_cpu = frequency;
    cudaMalloc(&_frequency_gpu, sizeof(*frequency));

    // assume that the CollisionFrequency model doesn't contain pointers else where
    cudaMemcpy(_frequency_gpu, _frequency_cpu, sizeof(*frequency), cudaMemcpyHostToDevice);
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

        // compute thermal energy for this cell
        for (int vi = 0; vi < nv; vi++) {
            double v = min_v + (vi + 0.5) * dv;
            int v_index = vi*(nc+2);

            thermal_energy += (v - v_avg) * (v - v_avg) * phi[v_index + ci];
        }

        // macroscopic properties for this cell
        double density = mass * n_particles / volume;
        double temp = 2 * thermal_energy / (density * R);

        // collision frequency for this cell
        double mu = frequency->collision_frequency(temp, n_particles);

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

void BGK::collide(Field<double> &phi, Field<double> &residual, Domain &domain, Equation &equation) {
    (void) phi;
    (void) domain;
    (void) equation;
}
