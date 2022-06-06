#pragma once
#include <complex>
#include <vector>

//#include "src/domain/rect_domain.h"

#include "src/domain/rect_domain_copy.h"
#include "src/initial_condition/initial_condition.h"
#include "src/serial_solver/base_serial_solver.h"

class FERectSolver:public BaseSolver
{
public:
    FERectSolver() = default;    
    FERectSolver(
        std::function<double(double, double)> potential, 
        double g, 
        RectangularDomain *domain);
    void solve();
    void solve_single_time(int k);
    double get_potential_value(int i, int j);

protected:
    RectangularDomain * domain;
    RectangularSpatialGrid potential_grid;
    void generate_potential_grid();
    std::complex<double> temporal_equation(int i, int j, int k);
};