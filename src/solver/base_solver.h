/**
 * @file base_solver.h
 * @author Minyoung Kim, Gyeonghun Kim
 * @brief Implementation file for base solver
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

#include "../domain/base_domain.h"
#include "../initial_condition/initial_condition.h"
/**
 * @brief Base Solver for Gross Piteavskill finite difference solver
 * 
 */
class BaseSolver
{
public:
    BaseSolver() = default;
    BaseSolver(float g);
    ~BaseSolver();

protected:
    float g;
    std::complex<float> temporal_equation(int i, int j, int k);
    std::string string_info;
};