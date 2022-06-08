#include "../../src/solver/serial_solver/crank_nicolson/cn_rect_solver.h"
#include "../../src/potential/harmonic_potential.h"
#include "../../src/utils.h"
#include <functional>
#include <iostream>
#include <complex>
#include "gtest/gtest.h"

TEST(CNSolverTest, InitializeSolveTest)
{
    bool all_passed = true;

    RectangularDomain *domain = new RectangularDomain(16, 16, 0, 1e-5, 2, -10, 10, -10, 10);

    auto initial_cond_function = [](float x, float y)
    { return std::complex<float>{1 * std::exp(-(x * x + y * y) / (9))}; };
    auto *initial_condition = new InitialCondition(initial_cond_function);
    initial_condition->assign_to_domain(domain);
    auto *potential = new HarmonicPotential(3, 5);
    potential->calcualte_potential_in_grid(domain);
    float g = -1;
    CNRectSolver solver = CNRectSolver(g, domain);
    
    solver.solve(1e-11, 101, "", false, false);
    ASSERT_TRUE(all_passed);
}
