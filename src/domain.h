#ifndef __DOMAIN_H_
#define __DOMAIN_H_

#include <json.hpp>
using json = nlohmann::json;

const unsigned BLOCK_SIZE = 256;


struct Domain {
public:
    Domain() {}
    Domain(double length, unsigned number_cells)
        : _length(length), _number_cells(number_cells) {}

    Domain(json json_data) {
        _length = json_data.at("length");
        _number_cells = json_data.at("number_cells");
    }

    double length() const {return _length;}
    double number_cells() const {return _number_cells;}
    double dx() const {return _length / _number_cells;}

    unsigned block_size() const {
        return BLOCK_SIZE;
    }

    unsigned number_blocks() const {
        return (_number_cells + BLOCK_SIZE -1) / BLOCK_SIZE;
    }

private:
    double _length;
    unsigned _number_cells;
};

#endif
