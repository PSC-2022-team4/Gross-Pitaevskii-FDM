#include "base_potential.h"

void BasePotential::calcualte_potential_in_grid(RectangularDomain *domain)
{
    int num_grid_1 = domain->get_num_grid_1();
    int num_grid_2 = domain->get_num_grid_2();
    for (auto i = 0; i < num_grid_1; ++i)
    {
        for (auto j = 0; j < num_grid_2; ++j)
        {
            auto point = domain->potential_grid->at(i, j);
            //Assign potential value in potential grid 
            point->value = {this->potential_function(point->x, point->y)};
        }
    }
}

std::string BasePotential::get_name()
{
    return this->name;
}

float BasePotential::potential_function(float x, float y)
{
    return 0;
}