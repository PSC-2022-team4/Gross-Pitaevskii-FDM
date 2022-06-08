#pragma once
#include <complex>
#include <vector>
#include <functional>
#include "src/domain/base_domain.h"
#include "src/initial_condition/initial_condition.h"
class BaseSolver
{
public:
    BaseSolver() = default;
    BaseSolver(float g);

protected:
    //  std::function<float(float, float)> potential_func;
    float g;
    std::complex<float> temporal_equation(int i, int j, int k);
    std::string string_info;
};