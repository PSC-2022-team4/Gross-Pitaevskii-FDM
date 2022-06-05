#pragma once
#include <complex>
#include <vector>
#include <functional>
#include <iostream>
#include <cmath>
#include <string.h>

#include "src/domain/rect_domain.h"
#include "src/initial_condition/initial_condition.h"
#include "src/serial_solver/base_serial_solver.h"
#include "src/serial_solver/forward_euler/fe_rect_solver.h"
class CNRectSolver : BaseSolver
{
public:
    CNRectSolver() = default;
    CNRectSolver(
        std::function<double(double, double)> potential, 
        double g, 
        RectangularDomain *domain
    );
    void generateRectangularDomain();
    void solve(double tolerance, int max_iter);
    void solve_single_time(int k, double tolerance, int max_iter);

protected:
    RectangularDomain *domain;
    RectangularSpatialGrid potential_grid;
    RectangularSpatialGrid *guess;
    FERectSolver *fe_solver;
    void generate_potential_grid();
    std::complex<double> temporal_equation(int i, int j, int k);
    std::complex<double> temporal_equation_from_guess(int i, int j);
    void initialize_guess_with_forward_euler(int k);
    void update_guess(int i, int j, int k);
    double calculate_error(int k);
};