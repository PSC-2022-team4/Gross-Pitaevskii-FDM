#pragma once
#include <complex>
#include <vector>
#include <functional>
#include "../../../domain/rect_domain.h"
#include "../../../initial_condition/initial_condition.h"
#include "../../base_solver.h"

class FERectSolver : public BaseSolver
{
public:
    FERectSolver() = default;
    FERectSolver( // std::function<float(float, float)> potential,
        float g,
        RectangularDomain *domain);
    void solve();
    void solve_single_time(int k);

protected:
    RectangularDomain *domain;
    // RectangularSpatialGrid potential_grid;
    void generate_potential_grid();
    float get_potential_value(int i, int j);
    std::complex<float> temporal_equation(int i, int j, int k);
};