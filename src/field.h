#ifndef FIELD_H_
#define FIELD_H_

#include <assert.h>
#include <stdexcept>
#include <iostream>

template <class T>
struct Field {
    Field () {}

    Field (unsigned length, unsigned nghost, unsigned components, bool gpu)
        : _length(length), _nghost(nghost), _components(components), _gpu(gpu)
    {
        if (_gpu) {
            auto code = cudaMalloc(&_data, (_length+2*_nghost) * _components * sizeof(T)); 

            if (code != cudaSuccess) {
                std::cerr << "cudaMalloc failure: " << cudaGetErrorString(code) << std::endl;
                throw new std::runtime_error("");
            }
        } 
        else {
            _data = new T[(_length+2*_nghost)*_components];
        }
    }

    Field (unsigned length, unsigned components, bool gpu) : Field(length, 0, components, gpu) {}
    Field (unsigned length, unsigned nghost, unsigned components) : Field(length, nghost, components, true) {}
    Field (unsigned length, unsigned components) : Field(length, 0, components, true) {}
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
        return _data[i + _nghost];
    }

    T& operator() (unsigned i, unsigned comp) {
        assert(i < _length && comp < _components);
        return _data[(_length+2*_nghost) * comp + i + _nghost];
    }

    T& flat_index(unsigned i) {
        return _data[i];
    }

    unsigned length () const { return _length; }
    unsigned number_components () const { return _components; }
    unsigned size() const { return _components * (_length+2*_nghost); }
    std::size_t memory() const { return size() * sizeof(T); }
    T* data() { return _data; }

private:
    T * _data = nullptr;
    unsigned _length;
    unsigned _nghost;
    unsigned _components;
    bool _gpu;
};


#endif
