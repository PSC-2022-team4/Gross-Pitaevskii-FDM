#include "../../src/solver/base_solver.h"
#include "../../src/solver/serial_solver/forward_euler/fe_rect_solver.h"

#include "../../src/potential/harmonic_potential.h"
#include "../../src/utils.h"
#include <iostream>
#include <complex>
#include "gtest/gtest.h"

TEST(FESolverTest, InitializeTest)
{
    // std::function<float(float, float)> potential;
    float g;
    RectangularDomain *domain = (new RectangularDomain(21, 21, 0, 1e-3, 3, -5, 5, -5, 5));
    auto initial_cond_function = [](float x, float y)
    { return std::complex<float>{1 * std::exp(-(x * x + y * y) / (9))}; };

    auto *initial_condition = new InitialCondition(initial_cond_function);
    initial_condition->assign_to_domain(domain);

    auto *potential = new HarmonicPotential(3, 5);
    potential->calcualte_potential_in_grid(domain);

    g = 1.;

    FERectSolver solver = FERectSolver(g, domain);

    ASSERT_FLOAT_EQ((*domain).at(10, 10, 0)->x, 0.);
    ASSERT_FLOAT_EQ((*domain).at(10, 10, 0)->x, 0.);
    ASSERT_FLOAT_EQ((*domain).at(10, 10, 0)->value.real(), 0.26607803);
    ASSERT_FLOAT_EQ((*domain).at(10, 10, 0)->value.imag(), 0.);
}
TEST(FESolverTest, SolveTest)
{
    float g;
    RectangularDomain *domain = (new RectangularDomain(21, 21, 0, 1e-3, 2, -5, 5, -5, 5));
    auto initial_cond_function = [](float x, float y)
    { return std::complex<float>{1 * std::exp(-(x * x + y * y) / (9))}; };

    auto *initial_condition = new InitialCondition(initial_cond_function);
    initial_condition->assign_to_domain(domain);

    auto *potential = new HarmonicPotential(3, 5);
    potential->calcualte_potential_in_grid(domain);

    g = 1.;

    FERectSolver solver = FERectSolver(g, domain);
    solver.solve("", false, false);
    ASSERT_FLOAT_EQ((*domain).at(10, 10, 1)->value.real(), 0.2660563);
    ASSERT_FLOAT_EQ((*domain).at(10, 10, 1)->value.imag(), -0.00013545621);
}