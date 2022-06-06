#pragma once
#include <complex>
#include <vector>
#include <functional>
#include "src/initial_condition/initial_condition.h"
class BaseSolver
{
public:
    BaseSolver() = default;
    BaseSolver(std::function<double(double, double)> potential, double g);

protected:
    std::function<double(double, double)> potential_func;
    double g;
    std::complex<double> temporal_equation(int i, int j, int k);
    std::string string_info;
};