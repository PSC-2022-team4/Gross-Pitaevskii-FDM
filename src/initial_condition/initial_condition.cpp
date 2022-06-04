#include "initial_condition.h"
#include <iostream>
#include <complex>
InitialCondition::InitialCondition(std::function<std::complex<double>(double, double)> initial_condition_function)
{
    this->initial_condition_function = initial_condition_function;
}

void InitialCondition::assign_to_domain(BaseDomain *domain)
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