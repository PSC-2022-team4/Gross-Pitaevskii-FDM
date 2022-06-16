/**
 * @file base_potential.h
 * @author Gyeonghun Kim, Minyoung Kim
 * @brief Header for the base potential class. Heritate this class for create custom potential shape. 
 * @version 0.1
 * @date 2022-06-08
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#pragma once

#include <iostream>
#include <string>

#include "../domain/rect_domain.h"

class BasePotential
{
public:
    BasePotential() = default;
    void calcualte_potential_in_grid(RectangularDomain *domain);
    std::string get_name();
    float potential_function(float, float);

protected:
    std::string name;
};