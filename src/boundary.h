#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

#include "field.h"

enum class BoundaryType {
    Periodic,
    Dirichlet,
    Neumann,
};

struct Domain;
struct Equation;

void fill_boundaries(Field<double> &phi, Domain &domain, Equation &equation);

BoundaryType boundary_type_from_string(std::string type);
std::string string_from_boundary_type(BoundaryType type);

#endif
