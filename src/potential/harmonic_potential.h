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
