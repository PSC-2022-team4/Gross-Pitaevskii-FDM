#include "src/serial_solver/base_serial_solver.h"
#include "src/serial_solver/forward_euler/fe_rect_solver.h"
#include "src/utils.h"
#include <iostream>
#include <complex>


bool test_fe_rect_solver(){
    bool all_passed = true;
    std::function<double(double, double)> potential;
    double g;
    RectangularDomain *domain = (new RectangularDomain(1001, 1001, 0, 10, 1001, -5, 5, -5, 5));
    auto initial_cond_function = [](double x, double y)
    { return std::complex<double>{1*std::exp(-(x*x + y*y)/(9))}; };

    auto *initial_condition =new  InitialCondition(initial_cond_function);
    initial_condition-> assign_to_domain(domain);
    
    
    potential= [](double x, double y ){
        return (double) 0.5 * (x*x + y *y);  };
    g = 1.; 
    
    
    FERectSolver solver = FERectSolver(potential, g, domain);
    
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
    if (!is_close((*domain).at(10, 10, 0)->wave_function.imag(), 0., 1e-12))
    {
        all_passed = false;
    }
    
    solver.solve();
    //At first step, psi(0,0) = 1+ i (4 e^-1 - 5) 
    double real = 1. ; 
    double imag = 4 * std::exp(-1./9.) - 5; 
    if (!is_close((*domain).at(10, 10, 1)->wave_function.real(), real, 1e-12))
    {
        all_passed = false;
    }if (!is_close((*domain).at(10, 10, 1)->wave_function.imag(), imag, 1e-12))
    {
        all_passed = false;
    }

    return all_passed;
    
}

bool test_all_fe_rect_solver()
{
    if (test_fe_rect_solver())
    {
        std::cout << "test_all_forward_euler_rectangular_serial_solver succeeded!" << std::endl;
    }
    else
    {
        std::cout << "test_all_forward_euler_rectangular_serial_solver failed!" << std::endl;
    }
}