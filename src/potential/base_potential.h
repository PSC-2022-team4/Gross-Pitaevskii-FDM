#pragma once
#include "../domain/rect_domain.h"
#include <string>
class BasePotential
{
public:
    BasePotential() = default;
    void calcualte_potential_in_grid(RectangularDomain *domain);
    std::string get_name();
    float potential_function(float, float);

protected:
    std::string name;
};