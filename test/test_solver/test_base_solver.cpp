#include "../../src/solver/base_solver.h"
#include "../../src/potential/harmonic_potential.h"
#include "../../src/utils.h"
#include <functional>
#include <iostream>
#include <complex>

bool test_base_solver()
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

    return all_passed;
}

bool test_all_base_solver()
{
    if (test_base_solver())
    {
        std::cout << "test_base_solver succeeded!" << std::endl;
    }
    else
    {
        std::cout << "test_base_solver failed!" << std::endl;
    }
}