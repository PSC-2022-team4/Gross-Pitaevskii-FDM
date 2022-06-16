/**
 * @file base_potential.cpp
 * @author Gyeonghun Kim, Minyoung Kim
 * @brief Implementation of the methods in potential class.
 * @version 0.1
 * @date 2022-06-08
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include "base_potential.h"

/**
 * @brief Calculate potential values in each grid based on domain point coordinate and potential function.
 * 
 * @param domain 
 */
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
            point->value = std::complex<float>{this->potential_function(point->x, point->y), 0};
        }
    }
    
}

/**
 * @brief Getter for name class variable
 * 
 * @return std::string 
 */
std::string BasePotential::get_name()
{
    return this->name;
}

/**
 * @brief Default potental function for constant zero potential.
 * 
 * @param x 
 * @param y 
 * @return float 
 */
float BasePotential::potential_function(float x, float y)
{
    return 0;
}