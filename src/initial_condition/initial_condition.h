#pragma once
#include "src/domain/base_domain.h"
#include <functional>
#include <complex>
class InitialCondition
{
    private:
        std::function<std::complex<double>(double, double)> initial_condition_function;
    public:
        InitialCondition() = default;
        InitialCondition(std::function<std::complex<double>(double, double)> initial_condition_function);
        void assign_to_domain(BaseDomain *domain);
};
