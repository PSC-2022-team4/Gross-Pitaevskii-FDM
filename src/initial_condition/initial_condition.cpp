/**
 * @file initial_condition.cpp
 * @author Gyeonghun Kim, Minyoung Kim
 * @brief Implementation of initial condition class
 * @version 0.1
 * @date 2022-06-04
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include "initial_condition.h"

/**
 * @brief Construct a new Initial Condition:: Initial Condition object
 * 
 * @param initial_condition_function 
 */
InitialCondition::InitialCondition(std::function<std::complex<float>(float, float)> initial_condition_function)
{
    this->initial_condition_function = initial_condition_function;
}

/**
 * @brief Assign initial condition to domain
 * 
 * @param domain 
 */
void InitialCondition::assign_to_domain(RectangularDomain *domain)
{

    for (auto i = 0; i < domain->get_num_grid_1(); ++i)
    {
        for (auto j = 0; j < domain->get_num_grid_2(); ++j)
        {
            auto point_data = domain->at(i, j, 0);
            domain->assign_initial_value(i, j, this->initial_condition_function(point_data->x, point_data->y));
        }
    }
    domain->normalize(0);
}