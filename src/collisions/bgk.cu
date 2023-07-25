#include "operator.h"
#include <json.hpp>

BGK::BGK(CollisionFrequency *frequency) {
    _frequency_cpu = frequency;
    cudaMalloc(&_frequency_gpu, sizeof(*frequency));
}

BGK::~BGK() {
    delete _frequency_cpu;
    cudaFree(_frequency_gpu);
}

void BGK::collide(Field<double> &phi, Domain &domain, Equation &equation) {
    (void) phi;
    (void) domain;
    (void) equation;
}
