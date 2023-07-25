#include "frequency.h"
#include "../constants.h"
#include <cmath>


__host__ __device__ 
HardSphere::HardSphere(double radius, double mass) 
    : _r_squared(radius * radius), _R(kB / mass) {}

__host__ __device__
HardSphere::HardSphere(const HardSphere &hard_sphere) {
    _r_squared = hard_sphere._r_squared;
    _R = hard_sphere._R;
}

__host__ __device__ 
double HardSphere::collision_frequency(double temp, double n) {
    return 16 * PI * n * _r_squared * sqrt(_R * temp);
}
