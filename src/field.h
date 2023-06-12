#ifndef __FIELD_H_
#define __FIELD_H_

#include <assert.h>
#include <stdexcept>
#include <iostream>

template <class T>
struct Field {
    Field () {}

    Field (unsigned length, unsigned components, bool gpu)
        : _length(length), _components(components), _gpu(gpu)
    {
        if (_gpu) {
            auto code = cudaMalloc(&_data, _length*_components*sizeof(T)); 

            if (code != cudaSuccess) {
                std::cerr << "cudaMalloc failure: " << cudaGetErrorString(code) << std::endl;
            }
        } 
        else {
            _data = new T[_length*_components];
        }
    }

    Field (unsigned length, unsigned components) : Field(length, components, true) {}
    Field (unsigned length, bool gpu) : Field(length, 1, gpu) {}
    Field (unsigned length) : Field(length, 1, true) {}

    ~Field() {
        if (_gpu) {
            cudaFree(_data);
        }
        else {
            delete [] _data;
        }
    }

    T& operator() (unsigned i) {
        assert(i < _length && _components == 1);
        return _data[i];
    }

    T& operator() (unsigned i, unsigned j) {
        assert(i < _length && j < _components);
        return _data[i*_length + j];
    }

    unsigned length () const { return _length; }
    unsigned number_components () const { return _components; }
    unsigned size() const { return _components * _length; }
    std::size_t memory() const { return size() * sizeof(T); }
    T* data() { return _data; }

private:
    T * _data;
    unsigned _length;
    unsigned _components;
    bool _gpu;
};

#endif
