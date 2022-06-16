#include "../../src/solver/base_solver.h"
#include "../../src/solver/parallel_solver/forward_euler/fe_rect_psolver.cuh"

#include "../../src/potential/harmonic_potential.h"
#include "../../src/utils.h"
#include <iostream>
#include <complex>

#include "gtest/gtest.h"

TEST(FEPSolverTest, InitializeSolveTest)
{
    bool all_passed = true;
    // std::function<float(float, float)> potential;

    float g;
    RectangularDomain *domain = (new RectangularDomain(32, 32, 0, 1e-4, 2, -5, 5, -5, 5));
    auto initial_cond_function = [](float x, float y)
    { return std::complex<float>{1e-10}; };

    auto *initial_condition = new InitialCondition(initial_cond_function);
    initial_condition->assign_to_domain(domain);

    auto *potential = new HarmonicPotential(3, 5);
    potential->calcualte_potential_in_grid(domain);
    g = -1.;
    FERectPSolver solver = FERectPSolver(g, domain, 0);

    solver.solve("", false, false);

    ASSERT_FLOAT_EQ((*domain).at(10, 10, 1)->value.real(), -0.027037166);
}
