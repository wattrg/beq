#ifndef _CONFIG_H
#define _CONFIG_H

#include <string>

struct Config {
    Config(std::string);

    int number_cells;
    double length;
    double velocity;
    int max_step;
    double max_time;
    int print_frequency;
    int plot_frequency;
    double cfl;
};

#endif
