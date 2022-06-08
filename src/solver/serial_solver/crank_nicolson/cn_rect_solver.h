#pragma once
#include <complex>
#include <vector>
#include <functional>
#include <iostream>
#include <cmath>
#include <string.h>

#include "../../../domain/rect_domain.h"
#include "../../../initial_condition/initial_condition.h"
#include "../../base_solver.h"
#include "../forward_euler/fe_rect_solver.h"
class CNRectSolver : public BaseSolver
{
public:
    CNRectSolver() = default;
    CNRectSolver(
        // std::function<float(float, float)> potential,
        float g,
        RectangularDomain *domain);
    void generateRectangularDomain();
    void solve(float tolerance, int max_iter);
    void solve_single_time(int k, float tolerance, int max_iter);

protected:
    RectangularDomain *domain;
    // RectangularSpatialGrid potential_grid;
    RectangularSpatialGrid *guess;
    FERectSolver *fe_solver;
    // void generate_potential_grid();
    std::complex<float> temporal_equation(int i, int j, int k);
    std::complex<float> temporal_equation_from_guess(int i, int j);
    void initialize_guess_with_forward_euler(int k);
    void update_guess(int i, int j, int k);
    float calculate_error(int k);
};