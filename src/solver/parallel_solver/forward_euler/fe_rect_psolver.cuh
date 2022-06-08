#pragma once
#include <complex>
#include <vector>
#include <functional>
#include "../../../domain/rect_domain.h"
#include "../../../initial_condition/initial_condition.h"
#include "../../base_solver.h"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#define nTx 16
#define nTy 16

class FERectPSolver : public BaseSolver
{
public:
    FERectPSolver() = default;
    FERectPSolver( // std::function<float(float, float)> potential,
        float g,
        RectangularDomain *domain);
    void solve(std::string dir_name = "");
    void solve_single_time(int k);

protected:
    RectangularDomain *domain;
    // RectangularSpatialGrid potential_grid;
    //void generate_potential_grid();
    float get_potential_value(int i, int j);
    std::complex<float> temporal_equation(int i, int j, int k);
};

__global__ void fe_rect_cusolver(float *psi_old_real,
                                 float *psi_old_imag,
                                 float *psi_new_real,
                                 float *psi_new_imag,
                                 float *potential,
                                 int n_x,
                                 int n_y,
                                 float g,
                                 float h_x,
                                 float h_y,
                                 float tau);