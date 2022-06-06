#include "src/parallel_solver/crank_nicolson/cn_rect_psolver.cuh"
#include "src/test/test_psolver/test_cn_rect_psolver.cuh"
#include "src/utils.h"
#include <functional>
#include <iostream>
#include <complex>

bool test_cn_psolver()
{
    bool all_passed = true;

    RectangularDomain *domain = (new RectangularDomain(101, 101, 0, 5, 1000, -5, 5, -5, 5));

    auto initial_cond_function = [](float x, float y)
    { return std::complex<float>{1. * exp(-((x-2.5)*(x-2.5)+y*y)/(1)) }; };
    auto *initial_condition = new InitialCondition(initial_cond_function);
    initial_condition->assign_to_domain(domain);

    auto potential = [](float x, float y)
    {
        return (float)0.5 * (x * x + y * y);
    };
    float g = 1;
    CNRectPSolver solver = CNRectPSolver(potential, g, domain);

<<<<<<< HEAD
    solver.solve(1e-11, 101);
=======

    solver.solve(1e-13, 101);
>>>>>>> 2de3edd5ca6d644807779450d5ab49958aa59459

    return all_passed;
}

bool test_all_cn_rect_psolver()
{
    bool passed;
    test_normalize(&passed);
    if (passed)
    {
        std::cout << "test_normalize succeeded!" << std::endl;
    }
    else
    {
        std::cout << "test_normalize failed!" << std::endl;
    }

    test_error_calculation(&passed);
    if (passed)
    {
        std::cout << "test_error_calculation succeeded!" << std::endl;
    }
    else
    {
        std::cout << "test_error_calculation failed!" << std::endl;
    }
    if (test_cn_psolver())
    {
        std::cout << "test_crank_nicolson_solver_creation succeeded!" << std::endl;
    }
    else
    {
        std::cout << "test_crank_nicolson_solver_creation failed!" << std::endl;
    }
}