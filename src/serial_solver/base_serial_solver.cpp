#include "base_serial_solver.h"

/**
 * @brief Construct a new Base Solver:: Base Solver object
 * 
 * @param potential_func_ 
 * @param g_ 
 */
BaseSolver::BaseSolver(std::function<double(double, double)> potential_func_, double g_)
    : potential_func(potential_func_), g(g_){};

/**
 * @brief defualt temporal equation 
 * 
 * @return std::complex<double> 0+0i
 */
std::complex<double> BaseSolver::temporal_equation(int i, int j, int k){
    return std::complex<double> {0};
}