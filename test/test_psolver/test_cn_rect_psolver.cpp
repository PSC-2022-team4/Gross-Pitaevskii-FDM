#include "src/solver/parallel_solver/crank_nicolson/cn_rect_psolver.cuh"
#include "src/test/test_psolver/test_cn_rect_psolver.cuh"
#include "src/potential/harmonic_potential.h"
#include "src/utils.h"
#include "src/potential/harmonic_potential.h"
#include <functional>
#include <iostream>
#include <complex>

bool test_cn_psolver()
{
    bool all_passed = true;

    RectangularDomain *domain = (new RectangularDomain(256, 256, 0, 5, 1000, -10, 10, -10, 10));

    auto initial_cond_function = [](float x, float y)
    { return std::complex<float>{1. * expf(-((x)*(x)+y*y)/(1)) }; };
    auto *initial_condition = new InitialCondition(initial_cond_function);
    initial_condition->assign_to_domain(domain);

    auto * potential = new HarmonicPotential(3, 5);
    potential->calcualte_potential_in_grid(domain);

    // auto potential = [](float x, float y)
    // {
    //     return 0.5 * ((x - 3) * (x - 3) * (x + 3) * (x + 3) + (y - 3) * (y - 3) * (y + 3) * (y + 3));
    // };
    float g = 1;
    CNRectPSolver solver = CNRectPSolver(g, domain);

    solver.solve(1e-11, 101);


    return all_passed;
}

bool test_all_cn_rect_psolver()
{
    bool passed;
    // test_normalize(&passed);
    // if (passed)
    // {
    //     std::cout << "test_normalize succeeded!" << std::endl;
    // }
    // else
    // {
    //     std::cout << "test_normalize failed!" << std::endl;
    // }

    // test_error_calculation(&passed);
    // if (passed)
    // {
    //     std::cout << "test_error_calculation succeeded!" << std::endl;
    // }
    // else
    // {
    //     std::cout << "test_error_calculation failed!" << std::endl;
    // }
    if (test_cn_psolver())
    {
        std::cout << "test_crank_nicolson_solver_creation succeeded!" << std::endl;
    }
    else
    {
        std::cout << "test_crank_nicolson_solver_creation failed!" << std::endl;
    }
}