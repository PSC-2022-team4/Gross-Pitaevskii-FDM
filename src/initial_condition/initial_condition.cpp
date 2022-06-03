#include "initial_condition.h"

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
            auto x = domain->get_domain_data()[0].at(i, j).x;
            auto y = domain->get_domain_data()[0].at(i, j).y;
            domain->get_domain_data()[0].at(i, j).wave_function = this->initial_condition_function(x, y);
        }
    }
}