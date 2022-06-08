#include "../../src/solver/serial_solver/crank_nicolson/cn_rect_solver.h"
#include "../../src/potential/harmonic_potential.h"
#include "../../src/utils.h"
#include <functional>
#include <iostream>
#include <complex>

bool test_cn_solver_creation()
{
    bool all_passed = true;

    RectangularDomain *domain = new RectangularDomain(256, 256, 0, 5, 1000, -10, 10, -10, 10);

    auto initial_cond_function = [](float x, float y)
    { return std::complex<float>{1 * std::exp(-(x * x + y * y) / (9))}; };
    auto *initial_condition = new InitialCondition(initial_cond_function);
    initial_condition->assign_to_domain(domain);
    auto *potential = new HarmonicPotential(3, 5);
    potential->calcualte_potential_in_grid(domain);
    float g = -1;
    CNRectSolver solver = CNRectSolver(g, domain);
    solver.solve(1e-11, 101);

    return all_passed;
}

bool test_all_cn_rect_solver()
{
    if (test_cn_solver_creation())
    {
        std::cout << "test_crank_nicolson_solver_creation succeeded!" << std::endl;
    }
    else
    {
        std::cout << "test_crank_nicolson_solver_creation failed!" << std::endl;
    }
}
