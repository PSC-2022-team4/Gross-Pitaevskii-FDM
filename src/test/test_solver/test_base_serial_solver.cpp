#include "src/serial_solver/base_serial_solver.h"
#include "src/utils.h"
#include <functional>
#include <iostream>
#include <complex>

bool test_base_serial_solver()
{
    bool all_passed = true;
    auto initial_cond_function = [](float x, float y)
    { return std::complex<float>{1 * std::exp(-(x * x + y * y) / (9)), 0}; };

    auto initial_condition = new InitialCondition(initial_cond_function);

    std::function<float(float, float)> potential = [](float x, float y)
    {
        return (float)0.5 * (x * x + y * y);
    };
    float g = 1;
    BaseDomain *domain = new BaseDomain(11, 11, 0, 10, 11);
    try{
        BaseSolver solver = BaseSolver(potential, g);
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