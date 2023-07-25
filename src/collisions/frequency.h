#ifndef COLLISION_FREQUENCY_H_
#define COLLISION_FREQUENCY_H_


#include <json.hpp>

using json = nlohmann::json;

class CollisionFrequency {
public:
    __host__ __device__ virtual ~CollisionFrequency () {};
    __host__ __device__ CollisionFrequency(){};

    __host__ __device__ virtual double collision_frequency(double temp, double n) = 0;
};

class HardSphere : public CollisionFrequency {
public:
    __host__ __device__ ~HardSphere() {};
    __host__ __device__ HardSphere(const HardSphere &hard_sphere);
    __host__ __device__ HardSphere(double radius, double mass);
    __host__ __device__ double collision_frequency(double temp, double n);

private:
    double _r_squared;
    double _R;
};

CollisionFrequency * make_collision_frequency(json frequency, json gas_model);


#endif
