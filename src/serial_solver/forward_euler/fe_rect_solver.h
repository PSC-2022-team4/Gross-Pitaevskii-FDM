#pragma once
#include <complex>
#include <vector>

//#include "src/domain/rect_domain.h"

#include "src/domain/rect_domain.h"
#include "src/initial_condition/initial_condition.h"
#include "src/serial_solver/base_serial_solver.h"

class FERectSolver:public BaseSolver
{
public:
    FERectSolver() = default;    
    FERectSolver(
        std::function<float(float, float)> potential, 
        float g, 
        RectangularDomain *domain);
    void solve();
    void solve_single_time(int k);
    float get_potential_value(int i, int j);
    void update_time();

protected:
    RectangularDomain * domain;
    RectangularSpatialGrid potential_grid;
    void generate_potential_grid();
    std::complex<float> temporal_equation(int i, int j, int k);
};