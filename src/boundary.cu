#include "boundary.h"
#include "equations/equation.h"

__global__
void fill_dirichlet(double *phi, double *values, int ci, int nc, int nv) {
    if (threadIdx.x == 0){
        for (int vi = 0; vi < nv; vi++) {
            phi[vi*nc + ci] = values[vi];
        }
    }
}

__global__
void copy_cell(double *phi, int src_cell, int dest_cell, int nc, int nv) {
    // phi: the data buffer
    // src_cell: the cell to copy
    // dest_cell: the cell to copy the data to
    // nc: the number of cells
    // nv: the number of velocities

    if (threadIdx.x == 0) {
        for (int vi = 0; vi < nv; vi++) {
            phi[vi*nc + dest_cell] = phi[vi*nc + src_cell];
        }
    }
}

void fill_boundaries(Field<double> &phi, Domain &domain, Equation &equation) {
    const int nc = domain.number_cells();
    const int nv = equation.number_components();

    // left boundary
    switch (domain.left_boundary()) {
        case BoundaryType::Periodic:
            copy_cell<<<1,1>>>(phi.data(), nc, 0, nc, nv);
            break;
        case BoundaryType::Neumann:
            copy_cell<<<1,1>>>(phi.data(), 1, 0, nc, nv);
            break;
        case BoundaryType::Dirichlet:
            fill_dirichlet<<<1,1>>>(phi.data(), domain.left_boundary_value().data(), 0, nc, nv);
    }

    // right boundary
    switch (domain.left_boundary()) {
        case BoundaryType::Periodic:
            copy_cell<<<1,1>>>(phi.data(), 1, nc+1, nc, nv);
            break;
        case BoundaryType::Neumann:
            copy_cell<<<1,1>>>(phi.data(), nc, nc+1, nc, nv);
            break;
        case BoundaryType::Dirichlet:
            fill_dirichlet<<<1,1>>>(phi.data(), domain.left_boundary_value().data(), nc+1, nc, nv);
    }
}