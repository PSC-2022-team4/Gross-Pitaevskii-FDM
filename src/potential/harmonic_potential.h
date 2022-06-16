/**
 * @file harmonic_potential.h
 * @author Minyoung Kim, Gyeonghun Kim
 * @brief Header file for harmoni potential (V = 0.5 * (omega_x^2 * x^2 + omega_y^2 * y^2) 
 * @version 0.1
 * @date 2022-06-08
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#pragma once

#include <iostream>

#include "base_potential.h"



class HarmonicPotential : public BasePotential
{
public:
    HarmonicPotential(float omega_x, float omega_y);
    void calcualte_potential_in_grid(RectangularDomain *domain);
    float potential_function(float x, float y);

private:
    float omega_x, omega_y;
};
