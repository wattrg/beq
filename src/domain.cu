#include "domain.h"

Domain::~Domain() {
    if (_left_boundary == BoundaryType::Dirichlet) {
        delete _left_boundary_value;
    }

    if (_right_boundary == BoundaryType::Dirichlet) {
        delete _right_boundary_value;
    }
}

Domain::Domain(json json_data) {
    _length = json_data.at("length");
    _number_cells = json_data.at("number_cells");

    // read the boundary data
    _left_boundary = boundary_type_from_string(
        json_data.at("left_boundary").at("type")
    );
    
    if (_left_boundary == BoundaryType::Dirichlet) {
        std::vector<double> left_value = json_data.at("left_boundary").at("value");
        _left_boundary_value = new Field<double>(1, 0, left_value.size(), true);

        auto code = cudaMemcpy(
            _left_boundary_value->data(), left_value.data(), _left_boundary_value->memory(), cudaMemcpyHostToDevice
        );
        if (code != cudaSuccess) {
            std::cerr << "Cuda error initialising left boundary condition: " << cudaGetErrorString(code) << std::endl;
            throw new std::runtime_error("Encountered cuda error");
        }
    }

    _right_boundary = boundary_type_from_string(
        json_data.at("right_boundary").at("type")
    );

    if (_right_boundary == BoundaryType::Dirichlet) {
        std::vector<double> right_value = json_data.at("right_boundary").at("value");
        _right_boundary_value = new Field<double>(1, 0, right_value.size(), true);

        auto code = cudaMemcpy(
            _right_boundary_value->data(), right_value.data(), _right_boundary_value->memory(), cudaMemcpyHostToDevice
        );
        if (code != cudaSuccess) {
            std::cerr << "Cuda error initialising right boundary condition: " << cudaGetErrorString(code) << std::endl;
            throw new std::runtime_error("Encountered cuda error");
        }
    }
}
