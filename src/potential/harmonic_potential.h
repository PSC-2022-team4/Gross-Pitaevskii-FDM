#pragma once
#include "base_potential.h"

class HarmonicPotential : public BasePotential
{
public:
    HarmonicPotential(float omega_x, float omega_y);
    float potential_function(float x, float y);

private:
    float omega_x, omega_y;
};
