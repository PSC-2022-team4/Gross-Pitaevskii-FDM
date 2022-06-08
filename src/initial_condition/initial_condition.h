#pragma once
#include "../domain/rect_domain.h"
#include <functional>
#include <complex>
#include <iostream>
class InitialCondition
{
private:
    std::function<std::complex<float>(float, float)> initial_condition_function;

public:
    InitialCondition() = default;
    InitialCondition(std::function<std::complex<float>(float, float)> initial_condition_function);
    void assign_to_domain(RectangularDomain *domain);
};
