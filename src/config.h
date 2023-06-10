#ifndef _CONFIG_H
#define _CONFIG_H

#include <string>

struct Config {
    Config(std::string);

    int number_cells() const { return _number_cells; }
    double length() const { return _length; }
    double velocity() const { return _velocity; }
    int max_step() const { return _max_step; }
    double max_time() const { return _max_time; }
    int print_frequency() const { return _print_frequency; }
    int plot_frequency() const { return _plot_frequency; }
    double cfl() const { return _cfl; }

private:
    int _number_cells;
    double _length;
    double _velocity;
    int _max_step;
    double _max_time;
    int _print_frequency;
    int _plot_frequency;
    double _cfl;
};

#endif
