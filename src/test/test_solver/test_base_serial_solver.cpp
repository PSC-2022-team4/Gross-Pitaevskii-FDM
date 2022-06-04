#include "src/serial_solver/base_serial_solver.h"
#include "src/utils.h"
#include <functional>
#include <iostream>
#include <complex>

bool test_base_serial_solver()
{
    bool all_passed = true;
    auto initial_cond_function = [](double x, double y)
    { return std::complex<double>{1 * std::exp(-(x * x + y * y) / (9)), 0}; };

    auto initial_condition = new InitialCondition(initial_cond_function);

    std::function<double(double, double)> potential = [](double x, double y)
    {
        return (double)0.5 * (x * x + y * y);
    };
    double g = 1;
    BaseDomain *domain = new BaseDomain(11, 11, 0, 10, 11);
    try{
        BaseSolver solver = BaseSolver(initial_condition, potential, g);
    }
    catch (int expn){
        all_passed = false;
    }

    return all_passed;
}

bool test_all_base_serial_solver()
{
    if (test_base_serial_solver())
    {
        std::cout << "test_base_serial_solver succeeded!" << std::endl;
    }
    else
    {
        std::cout << "test_base_serial_solver failed!" << std::endl;
    }
}