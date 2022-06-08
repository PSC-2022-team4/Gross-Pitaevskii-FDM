#include "harmonic_potential.h"


HarmonicPotential::HarmonicPotential(float omega_x, float omega_y)
: BasePotential(), omega_x(omega_x), omega_y(omega_y)
{
    this->name = std::string("Harmonic");
}

float HarmonicPotential::potential_function(float x, float y)
{
    return 0.5 * (this->omega_x * x * x + this->omega_y * y * y);
}
