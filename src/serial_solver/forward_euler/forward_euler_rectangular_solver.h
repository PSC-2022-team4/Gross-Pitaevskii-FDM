#pragma once
#include <complex>
#include <vector>
#include <functional>
#include "src/domain/rectangular_domain.h"
#include "src/initial_condition/initial_condition.h"
#include "src/serial_solver/base_serial_solver.h"
class ForwardEulerRectangularSolver : BaseSolver
{
public:
    ForwardEulerRectangularSolver() = default;
    ForwardEulerRectangularSolver(InitialCondition initialCondition, std::function<double(double, double)> potential, double g, RectangularDomain* rectangularDomain);
    void applyInitialCondition();
    void solve();
    void solve_single_time(int k);

protected:
    InitialCondition initialCondition;
    RectangularDomain * domain;
    std::function<double(double, double)> potential_func;
    double g;
    std::complex<double> temporal_equation(int i, int j, int k);
    
};