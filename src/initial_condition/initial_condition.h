/**
 * @file initial_condition.h
 * @author Gyeonghun Kim, Minyoung Kim
 * @brief Header for the initial condition class
 * @version 0.1
 * @date 2022-06-04
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#pragma once
#include <functional>
#include <complex>
#include <iostream>
#include <iostream>
#include <complex>

#include "../domain/rect_domain.h"

/**
 * @brief Initial condition class
 * 
 */
class InitialCondition
{
private:
    std::function<std::complex<float>(float, float)> initial_condition_function;

public:
    InitialCondition() = default;
    InitialCondition(std::function<std::complex<float>(float, float)> initial_condition_function);
    void assign_to_domain(RectangularDomain *domain);
};
