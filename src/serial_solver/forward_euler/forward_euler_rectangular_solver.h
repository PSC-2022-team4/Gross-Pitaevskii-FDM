#pragma once
#include <complex>
#include <vector>
#include <functional>
#include "src/domain/base_domain.h"
#include "src/initial_condition/initial_condition.h"
#include "src/serial_solver/base_serial_solver.h"
class ForwardEulerRectangularSolver : BaseSolver
{
public:
    ForwardEulerRectangularSolver() = default;
    ForwardEulerRectangularSolver(InitialCondition initialCondition, std::function<double(double, double)> potential, double g);

protected:
    InitialCondition initialCondition;
    BaseDomain baseDomain;
    std::function<double(double, double)> potential_func;
    double g;
    std::complex<double> temporal_equation(double x, double y, double t);
    void update();
};