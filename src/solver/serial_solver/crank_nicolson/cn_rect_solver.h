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

class CNRectSolver : FERectSolver
{
public:
    CNRectSolver() = default;
    CNRectSolver(
        float g,
        RectangularDomain *domain);
    void solve(float tolerance, int max_iter, std::string dir_name = "");
    void solve_single_time(int k, float tolerance, int max_iter);

protected:
    RectangularSpatialGrid *guess;
    void initialize_guess_with_forward_euler(int k);
    void update_guess(int i, int j, int k);
    float calculate_error(int k);
};
