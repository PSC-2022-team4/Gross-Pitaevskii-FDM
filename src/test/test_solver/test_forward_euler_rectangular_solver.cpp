#include "src/serial_solver/forward_euler/forward_euler_rectangular_solver.h"
#include "src/utils.h"
#include <iostream>
#include <complex>


bool test_forward_euler_rectangular_solver(){
    bool all_passed = true;
    auto initial_cond_function = [](double x, double y)
    { return std::complex<double>{1*std::exp(-(x*x + y*y)/(9))}; };

    InitialCondition initial_condition = InitialCondition(initial_cond_function);
    
    std::function<double(double, double)> potential = [](double x, double y ){
        return (double) 0.5 * (x*x + y *y);  
    };
    double g = 1 ; 
    RectangularDomain* domain = new RectangularDomain(21, 21, 0, 10, 11, -10, 10, -10, 10);
    ForwardEulerRectangularSolver solver = ForwardEulerRectangularSolver(initial_condition, potential, g, domain);
    
    // if(!is_close((*domain).at(10,10, 0)->x , 0., 1e-12)){
    //     all_passed = false;
    // }
    // if (!is_close((*domain).at(10, 10, 0)->y, 0., 1e-12))
    // {
    //     all_passed = false;
    // }
    // if (!is_close((*domain).at(10, 10, 0)->wave_function.real(), 1., 1e-12))
    // {
    //     all_passed = false;
    // }
    // if (!is_close((*domain).at(10, 10, 0)->wave_function.imag(), 0., 1e-12))
    // {
    //     all_passed = false;
    // }
    // //At first step, psi(0,0) = 1+ i (4 e^-1 - 5) 
    // double real = 0. ; 
    // double imag = 4 * std::exp(-1) - 5; 
    // if (!is_close((*domain).at(10, 10, 1)->wave_function.real(), real, 1e-12))
    // {
    //     all_passed = false;
    // }if (!is_close((*domain).at(10, 10, 1)->wave_function.imag(), imag, 1e-12))
    // {
    //     all_passed = false;
    // }

    return all_passed;
    
}
