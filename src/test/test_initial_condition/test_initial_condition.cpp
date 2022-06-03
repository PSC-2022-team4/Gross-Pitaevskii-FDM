#include "src/domain/rectangular_domain.h"
#include "src/initial_condition/initial_condition.h"
#include "src/utils.h"
#include <iostream>
#include <complex>

bool test_initial_condition_rectangular()
{
    auto domain = RectangularDomain(101, 101, 0, 10, 11, -10, 10, -10, 10);
    auto initial_cond_function = [](double x, double y)
    { return std::complex<double>{1*std::exp(-(x*x + y*y)/(9))}; };

    auto initial_condition = InitialCondition(initial_cond_function);
    initial_condition.assign_to_domain(&domain);
    std::cout << "." << std::endl;
    std::cout << domain.get_domain_data()[0].at(51, 51).wave_function << std::endl;
    bool all_passed = true;

    // if (!is_close(domain.x, 0.1, 1e-12))
    // {
    //     all_passed = false;
    // }
    // if (!is_close(grid_point.y, 0.1, 1e-12))
    // {
    //     all_passed = false;
    // }
    // if (!is_close(grid_point.wave_function.real(), 10., 1e-12))
    // {
    //     all_passed = false;
    // }
    // if (!is_close(grid_point.wave_function.imag(), 1., 1e-12))
    // {
    //     all_passed = false;
    // }
    return all_passed;
}

bool test_all_initial_condition()
{
    if (test_initial_condition_rectangular())
    {
        std::cout << "test_initial_condition_rectangular succeeded!" << std::endl;
    }
    else
    {
        std::cout << "test_initial_condition_rectangular failed!" << std::endl;
    }
}