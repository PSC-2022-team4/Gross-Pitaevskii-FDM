#include "base_solver.h"

/**
 * @brief Construct a new Base Solver:: Base Solver object
 * 
 * @param g_ 
 */
BaseSolver::BaseSolver(float g_) 
    : g(g_){
        this -> string_info = "Base_Solver";
    };                   

BaseSolver::~BaseSolver(){};
/**
 * @brief defualt temporal equation 
 * 
 * @return std::complex<float> 0+0i
 */
std::complex<float> BaseSolver::temporal_equation(int i, int j, int k){
    return std::complex<float> {0};
}