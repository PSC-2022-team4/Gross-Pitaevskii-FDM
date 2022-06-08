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
    std::cout << "." << std::endl;
    RectangularDomain *domain = (new RectangularDomain(32, 32, 0, 1e-4, 2, -5, 5, -5, 5));
    std::cout << "." << std::endl;
    auto initial_cond_function = [](float x, float y)
    { return std::complex<float>{1e-10}; };

    auto *initial_condition = new InitialCondition(initial_cond_function);
    initial_condition->assign_to_domain(domain);

    auto *potential = new HarmonicPotential(3, 5);
    potential->calcualte_potential_in_grid(domain);
    g = -1.;
    FERectPSolver solver = FERectPSolver(g, domain, 0);

    solver.solve();

    ASSERT_FLOAT_EQ((*domain).at(10, 10, 1)->value.real(), 0.096875511);
    ASSERT_TRUE(abs((*domain).at(10, 10, 1)->value.imag()) < 1e-6);
}

bool test_fe_rect_psolver()
{
}

bool test_all_fe_rect_psolver()
{
    if (test_fe_rect_psolver())
    {
        std::cout << "test_all_forward_euler_rectangular_serial_solver succeeded!" << std::endl;
    }
    else
    {
        std::cout << "test_all_forward_euler_rectangular_serial_solver failed!" << std::endl;
    }
}