#ifndef __DOMAIN_H_
#define __DOMAIN_H_

#include <json.hpp>
#include <vector>

#include "boundary.h"

using json = nlohmann::json;


const unsigned BLOCK_SIZE = 256;

struct Domain {
public:
    Domain() {}
    // Domain(double length, unsigned number_cells)
    //     : _length(length), _number_cells(number_cells) {}

    Domain(json json_data) {
        _length = json_data.at("length");
        _number_cells = json_data.at("number_cells");

        // read the boundary data
        _left_boundary = boundary_type_from_string(
            json_data.at("left_boundary").at("type")
        );

        _right_boundary = boundary_type_from_string(
            json_data.at("right_boundary").at("type")
        );
    }

    double length() const {return _length;}
    double number_cells() const {return _number_cells;}
    double dx() const {return _length / _number_cells;}

    BoundaryType left_boundary() { return _left_boundary; }
    BoundaryType right_boundary() { return _right_boundary; }

    std::vector<double> left_boundary_value() { return _left_boundary_value; }
    std::vector<double> right_boundary_value() { return _right_boundary_value; }

    unsigned number_ghost () const { return _number_ghost; }

    unsigned block_size() const {
        return BLOCK_SIZE;
    }

    unsigned number_blocks() const {
        return (_number_cells + BLOCK_SIZE -1) / BLOCK_SIZE;
    }

private:
    double _length;
    unsigned _number_cells;
    unsigned _number_ghost = 1;

    BoundaryType _left_boundary;
    std::vector<double> _left_boundary_value;

    BoundaryType _right_boundary;
    std::vector<double> _right_boundary_value;
};

#endif
