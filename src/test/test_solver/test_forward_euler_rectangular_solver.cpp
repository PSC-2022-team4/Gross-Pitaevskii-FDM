#include "src/serial_solver/base_serial_solver.h"
#include "src/serial_solver/forward_euler/forward_euler_rectangular_solver.h"
#include "src/utils.h"
#include <iostream>
#include <complex>


bool test_forward_euler_rectangular_solver(){
    bool all_passed = true;
    InitialCondition initial_condition;
    std::function<double(double, double)> potential;
    double g; 
    RectangularDomain* domain;
    auto initial_cond_function = [](double x, double y)
    { return std::complex<double>{1*std::exp(-(x*x + y*y)/(9))}; };

    initial_condition = InitialCondition(initial_cond_function);
    
    potential= [](double x, double y ){
        return (double) 0.5 * (x*x + y *y);  
    };
    g = 1. ; 
    
    domain = (new RectangularDomain(21, 21, 0, 10, 11, -10, 10, -10, 10));
    
    try{
        ForwardEulerRectangularSolver solver = ForwardEulerRectangularSolver(initial_condition, potential, g, domain);
    }
    catch (const std::bad_function_call &e)
    {
        std::cout << e.what() << '\n';
        all_passed = false;
    }

    if(!is_close((*domain).at(10,10, 0)->x , 0., 1e-12)){
        all_passed = false;
    }
    if (!is_close((*domain).at(10, 10, 0)->y, 0., 1e-12))
    {
        all_passed = false;
    }
    if (!is_close((*domain).at(10, 10, 0)->wave_function.real(), 1., 1e-12))
    {
        all_passed = false;
    }
    if (!is_close((*domain).at(10, 10, 0)->wave_function.imag(), 0., 1e-12))
    {
        all_passed = false;
    }
    //At first step, psi(0,0) = 1+ i (4 e^-1 - 5) 
    double real = 0. ; 
    double imag = 4 * std::exp(-1) - 5; 
    if (!is_close((*domain).at(10, 10, 1)->wave_function.real(), real, 1e-12))
    {
        all_passed = false;
    }if (!is_close((*domain).at(10, 10, 1)->wave_function.imag(), imag, 1e-12))
    {
        all_passed = false;
    }

    return all_passed;
    
}

bool test_all_forward_euler_rectangular_serial_solver()
{
    if (test_forward_euler_rectangular_solver())
    {
        std::cout << "test_base_serial_solver succeeded!" << std::endl;
    }
    else
    {
        std::cout << "test_base_serial_solver failed!" << std::endl;
    }
}