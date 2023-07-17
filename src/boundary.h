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

#endif
