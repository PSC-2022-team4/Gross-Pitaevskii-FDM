/**
 * @file fe_rect_solver.h
 * @author Minyoung Kim, Gyeonghun Kim
 * @brief Header file for serial forward euler solver
 * @version 0.1
 * @date 2022-06-05
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#pragma once
#include <complex>
#include <vector>
#include <functional>
#include <iostream>
#include <string>

#include "../../../domain/rect_domain.h"
#include "../../../initial_condition/initial_condition.h"
#include "../../base_solver.h"

/**
 * @brief Forward Euler Serial Solver 
 * 
 */
class FERectSolver : public BaseSolver
{
public:
    FERectSolver() = default;    
    FERectSolver(               
        float g, 
        RectangularDomain *domain);
    void solve(std::string dir_name="", bool print_info=true, bool save_data=true);
    void solve_single_time(int k);
    //void update_time();

protected:
    RectangularDomain * domain;
    float get_potential_value(int i, int j);
    std::complex<float> temporal_equation(int i, int j, int k);
};