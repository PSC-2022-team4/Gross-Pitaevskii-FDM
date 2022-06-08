#include <iostream>
#include <mpi.h>
#include "domain/rect_domain.h"
#include "initial_condition/initial_condition.h"
#include "potential/harmonic_potential.h"
#include "solver/base_solver.h"
#include "solver/parallel_solver/forward_euler/fe_rect_psolver.cuh"
#include "solver/parallel_solver/crank_nicolson/cn_rect_psolver.cuh"
#include "solver/serial_solver/forward_euler/fe_rect_solver.h"
#include "solver/serial_solver/crank_nicolson/cn_rect_solver.h"

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    RectangularDomain *domain = (new RectangularDomain(256, 256, 0, 5, 1000, -10, 10, -10, 10));

    auto initial_cond_function = [](float x, float y)
    { return std::complex<float>{float(1.) * expf(-((x) * (x) + y * y) / (1))}; };
    auto *initial_condition = new InitialCondition(initial_cond_function);
    initial_condition->assign_to_domain(domain);

    auto *potential = new HarmonicPotential(3., 5.);
    potential->calcualte_potential_in_grid(domain);

    float g = 1;
    // CNRectSolver solver = CNRectSolver(g, domain);
    // FERectSolver solver = FERectSolver(g, domain);
    // CNRectPSolver solver = CNRectPSolver(g, domain, 0);
    FERectPSolver solver = FERectPSolver(g, domain, 0);
    
    // solver.solve(1e-11, 101);
    solver.solve();
    MPI_Finalize();
}
