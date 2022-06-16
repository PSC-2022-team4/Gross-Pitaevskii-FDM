/**
 * @file base_solver.cpp
 * @author Minyoung Kim, Gyeonghun Kim
 * @brief Header file for base solver
 * @version 0.1
 * @date 2022-06-05
 * 
 * @copyright Copyright (c) 2022
 * 
 */
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


/**
 * @brief Destroy the Base Solver:: Base Solver object
 * 
 */
BaseSolver::~BaseSolver(){};
/**
 * @brief defualt temporal equation 
 * 
 * @return std::complex<float> 0+0i
 */


std::complex<float> BaseSolver::temporal_equation(int i, int j, int k){
    return std::complex<float> {0};
}