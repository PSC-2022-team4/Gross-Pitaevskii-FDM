#include "../../src/solver/base_solver.h"
#include "../../src/potential/harmonic_potential.h"
#include "../../src/utils.h"
#include <functional>
#include <iostream>
#include <complex>
#include "gtest/gtest.h"

TEST(BaseSolverTest, InitializeTest)
{
    bool all_passed = true;

    BaseDomain *domain = new BaseDomain(11, 11, 0, 10, 11);

    auto initial_cond_function = [](float x, float y)
    { return std::complex<float>{1 * std::exp(-(x * x + y * y) / (9)), 0}; };
    auto initial_condition = new InitialCondition(initial_cond_function);
    auto *potential = new HarmonicPotential(3, 5);

    float g = 1;

    try
    {
        BaseSolver solver = BaseSolver(g);
    }
    catch (int expn)
    {
        all_passed = false;
    }

    ASSERT_TRUE(all_passed);
}