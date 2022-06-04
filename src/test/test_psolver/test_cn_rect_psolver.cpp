#include "src/parallel_solver/crank_nicolson/cn_rect_psolver.cuh"
#include "src/utils.h"
#include <functional>
#include <iostream>
#include <complex>

bool test_cn_psolver()
{
    bool all_passed = true;

    RectangularDomain *domain = (new RectangularDomain(41, 41, 0, 10, 1001, -5, 5, -5, 5));

    auto initial_cond_function = [](double x, double y)
    { return std::complex<double>{1 * std::exp(-(x * x + y * y) / (9))}; };
    auto *initial_condition = new InitialCondition(initial_cond_function);
    initial_condition->assign_to_domain(domain);

    auto potential = [](double x, double y)
    {
        return (double)0.5 * (x * x + y * y);
    };
    double g = -1;
    CNRectPSolver solver = CNRectPSolver(potential, g, domain);


    solver.solve(1e-13, 101);

    return all_passed;
}

bool test_all_cn_rect_psolver()
{
    if (test_cn_psolver())
    {
        std::cout << "test_crank_nicolson_solver_creation succeeded!" << std::endl;
    }
    else
    {
        std::cout << "test_crank_nicolson_solver_creation failed!" << std::endl;
    }
}
