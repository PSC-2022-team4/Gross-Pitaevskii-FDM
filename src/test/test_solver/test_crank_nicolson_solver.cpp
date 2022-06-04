#include "src/serial_solver/crank_nicolson/crank_nicolson_rectangular_solver.h"
#include "src/utils.h"
#include <functional>
#include <iostream>
#include <complex>

bool test_crank_nicolson_solver_creation()
{
    bool all_passed = true;
    auto initial_cond_function = [](double x, double y)
    { return std::complex<double>{1 * std::exp(-(x * x + y * y) / (9))}; };

    auto initial_condition = InitialCondition(initial_cond_function);

    std::function<double(double, double)> potential = [](double x, double y)
    {
        return (double)0.5 * (x * x + y * y);
    };
    double g = 1;
    RectangularDomain *domain = new RectangularDomain(21, 21, 0, 10, 11, -10, 10, -10, 10);
    ForwardEulerRectangularSolver solver = ForwardEulerRectangularSolver(initial_condition, potential, g, domain);

    try
    {
        BaseSolver solver = BaseSolver(initial_condition, potential, g, domain);
    }
    catch (int expn)
    {
        all_passed = false;
    }

    return all_passed;
}

bool test_all_crank_nicolson_solver()
{
    if (test_crank_nicolson_solver_creation())
    {
        std::cout << "test_crank_nicolson_solver_creation succeeded!" << std::endl;
    }
    else
    {
        std::cout << "test_crank_nicolson_solver_creation failed!" << std::endl;
    }
}