/**
 * @file harmonic_potential.cpp
 * @author Gyeonghun Kim, Minyoung Kim
 * @brief Header file of harmonic potential class
 * @version 0.1
 * @date 2022-06-08
 * 
 * @copyright Copyright (c) 2022
 * 
 */
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
void HarmonicPotential::calcualte_potential_in_grid(RectangularDomain *domain)
{
    int num_grid_1 = domain->get_num_grid_1();
    int num_grid_2 = domain->get_num_grid_2();
    for (auto i = 0; i < num_grid_1; ++i)
    {
        for (auto j = 0; j < num_grid_2; ++j)
        {
            auto point = domain->potential_grid->at(i, j);

            //Assign potential value in potential grid
            point->value = std::complex<float>{this->potential_function(point->x, point->y), 0};
        }
    }
}