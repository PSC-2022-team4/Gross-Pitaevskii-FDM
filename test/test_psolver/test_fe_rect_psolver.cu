#include "../../src/solver/base_solver.h"
#include "../../src/solver/parallel_solver/forward_euler/fe_rect_psolver.cuh"

#include "../../src/potential/harmonic_potential.h"
#include "../../src/utils.h"
#include <iostream>
#include <complex>

bool test_fe_rect_psolver()
{
    bool all_passed = true;
    // std::function<float(float, float)> potential;

    float g;
    std::cout << "." << std::endl;
    RectangularDomain *domain = (new RectangularDomain(101, 101, 0, 1e-2, 11, -5, 5, -5, 5));
    std::cout << "." << std::endl;
    auto initial_cond_function = [](float x, float y)
    { return std::complex<float>{1e-10}; };

    auto *initial_condition = new InitialCondition(initial_cond_function);
    initial_condition->assign_to_domain(domain);

    auto *potential = new HarmonicPotential(3, 5);
    potential->calcualte_potential_in_grid(domain);
    // potential= [](float x, float y ){
    //     return (float) 0.5 * 5 * (x*x + y *y);  };
    g = -1.;

    std::cout << "." << std::endl;
    FERectPSolver solver = FERectPSolver(g, domain, 0);

    if (!is_close((*domain).at(10, 10, 0)->x, 0., 1e-12))
    {
        all_passed = false;
    }
    if (!is_close((*domain).at(10, 10, 0)->y, 0., 1e-12))
    {
        all_passed = false;
    }
    if (!is_close((*domain).at(10, 10, 0)->value.real(), 1., 1e-12))
    {
        all_passed = false;
    }
    if (!is_close((*domain).at(10, 10, 0)->value.imag(), 0., 1e-12))
    {
        all_passed = false;
    }
    if (!is_close((*domain).at(10, 10, 0)->value.imag(), 0., 1e-12))
    {
        all_passed = false;
    }

    solver.solve();
    //At first step, psi(0,0) = 1+ i (4 e^-1 - 5)
    float real = 1.;
    float imag = 4 * std::exp(-1. / 9.) - 5;
    if (!is_close((*domain).at(10, 10, 1)->value.real(), real, 1e-12))
    {
        all_passed = false;
    }
    if (!is_close((*domain).at(10, 10, 1)->value.imag(), imag, 1e-12))
    {
        all_passed = false;
    }

    return all_passed;
}

bool test_all_fe_rect_psolver()
{
    if (test_fe_rect_psolver())
    {
        std::cout << "test_all_forward_euler_rectangular_serial_solver succeeded!" << std::endl;
    }
    else
    {
        std::cout << "test_all_forward_euler_rectangular_serial_solver failed!" << std::endl;
    }
}