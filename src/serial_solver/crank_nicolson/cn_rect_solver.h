#pragma once
#include <complex>
#include <vector>
#include <functional>
#include <iostream>
#include <cmath>
#include <string.h>
#include "src/domain/rect_domain.h"
#include "src/initial_condition/initial_condition.h"
#include "src/serial_solver/forward_euler/fe_rect_solver.h"

class CNRectSolver : FERectSolver
{
public:
    CNRectSolver() = default;
    CNRectSolver(
        std::function<float(float, float)> potential, 
        float g, 
        RectangularDomain *domain
    );
    void solve(float tolerance, int max_iter, std::string dir_name="");
    void solve_single_time(int k, float tolerance, int max_iter);

protected:
    RectangularSpatialGrid *guess;
    //FERectSolver *fe_solver;
    void initialize_guess_with_forward_euler(int k);
    void update_guess(int i, int j, int k);
    float calculate_error(int k);
};