#include "base_potential.h"

void BasePotential::calcualte_potential_in_grid(RectangularDomain *domain)
{
    int num_grid_1 = domain->get_num_grid_1();
    int num_grid_2 = domain->get_num_grid_2();
    float x_start = domain->get_x_start();
    float y_start = domain->get_y_start();
    float x_end = domain->get_x_end();
    float y_end = domain->get_y_end();;
    domain->potential_grid = RectangularSpatialGrid(num_grid_1, num_grid_2, x_start, x_end, y_start, y_end);
    for (auto i = 0; i < num_grid_1; ++i)
    {
        for (auto j = 0; j < num_grid_2; ++j)
        {
            auto point = domain->potential_grid.at(i, j);
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