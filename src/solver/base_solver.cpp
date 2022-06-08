#include "base_solver.h"

/**
 * @brief Construct a new Base Solver:: Base Solver object
 * 
 * @param potential_func_ 
 * @param g_ 
 */
BaseSolver::BaseSolver(float g_) // std::function<float(float, float)> potential_func_,
    : g(g_){};                   // potential_func(potential_func_),

/**
 * @brief defualt temporal equation 
 * 
 * @return std::complex<float> 0+0i
 */
std::complex<float> BaseSolver::temporal_equation(int i, int j, int k){
    return std::complex<float> {0};
}