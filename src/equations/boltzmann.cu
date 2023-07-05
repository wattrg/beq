#include "equation.h"

Boltzmann::Boltzmann(json json_data) {
    _min_v = json_data.at("min_v");
    _max_v = json_data.at("max_v");
    _nv = json_data.at("n_vel_increments");
    _dv = (_max_v - _min_v) / _nv;
}

__global__
void eval_boltzmann_residual(double *phi, double *residual, 
                             double dx, double dv, int nc, int nv, double min_v)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    for (int ci = index; ci < nc; ci += stride) {
        for (int vi = 0; vi < nv; vi++) {
            double phi_minus, phi_plus;
            double v = min_v + (vi + 0.5) * dv;

            int v_index = vi*nc;

            // left boundary
            if (ci == 0) {
                if (v > 0) {
                    // right moving particles
                    phi_plus = phi[v_index + ci];
                    phi_minus = phi[v_index + nc-1];
                }
                else if (v < 0){
                    // left moving particles
                    phi_plus = phi[v_index + ci + 1];
                    phi_minus = phi[v_index + ci];
                }
                else {
                    // stationary particles
                    phi_plus = phi[v_index + ci];
                    phi_minus = phi[v_index + ci];
                }
            }

            // right boundary
            if (ci == nc-1) {
                if (v < 0) {
                    // left moving particle
                    phi_plus = phi[v_index + 0];
                    phi_minus = phi[v_index + ci];
                }
                else if (v > 0) {
                    // right moving particle
                    phi_plus = phi[v_index + ci];
                    phi_minus = phi[v_index + ci - 1];
                }
                else {
                    // stationary particles
                    phi_plus = phi[v_index + ci];
                    phi_minus = phi[v_index + ci];
                }
            }

            // interior cell
            if (v < 0) {
                // left moving particles
                phi_plus = phi[v_index + ci + 1];
                phi_minus = phi[v_index + ci];
            }
            else if (v > 0) {
                // right moving particles
                phi_plus = phi[v_index + ci];
                phi_minus = phi[v_index + ci - 1];
            }
            else {
                // stationary particles
                phi_plus = phi[v_index + ci];
                phi_minus = phi[v_index + ci];
            }

            residual[vi*nc + ci] = v * (phi_plus - phi_minus);
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
    double max_v = fmax(fabs(_min_v), fabs(_max_v));
    return domain.dx() / max_v;
}
