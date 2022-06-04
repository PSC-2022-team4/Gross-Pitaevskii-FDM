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

    auto initial_condition = new InitialCondition(initial_cond_function);
    auto potential = [](double x, double y)
    {
        return (double)0.5 * (x * x + y * y);
    };
    double g = 1;
    RectangularDomain *domain = new RectangularDomain(21, 21, 0, 10, 11, -10, 10, -10, 10);

    CrankNicolsonRectangularSolver solver = CrankNicolsonRectangularSolver(potential, g, domain);

    solver.solve(1e-3, 101);

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

