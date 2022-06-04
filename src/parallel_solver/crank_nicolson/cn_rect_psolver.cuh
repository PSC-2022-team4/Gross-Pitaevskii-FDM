#pragma once
#include <complex>
#include <vector>
#include <functional>
#include <iostream>
#include <cmath>

#include "src/domain/rect_domain.h"
#include "src/initial_condition/initial_condition.h"
#include "src/serial_solver/base_serial_solver.h"
#include "src/serial_solver/forward_euler/fe_rect_solver.h"

#include <cuda.h>
#include <cuda_runtime.h>

#define nTx 16
#define nTy 16

class CNRectPSolver : BaseSolver
{
public:
    CNRectPSolver() = default;
    CNRectPSolver(
        std::function<double(double, double)> potential,
        double g,
        RectangularDomain *domain);
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

__global__ void cn_rect_cusolver(double *psi_old_real,
                                 double *psi_old_imag,
                                 double *psi_new_real_trial,
                                 double *psi_new_imag_trial,
                                 double *psi_new_real,
                                 double *psi_new_imag,
                                 double *potential,
                                 int n_x,
                                 int n_y,
                                 double g,
                                 double h_x,
                                 double h_y,
                                 double tau);

__global__ void calculate_error_on_device(double *psi_1_real,
                                          double *psi_1_imag,
                                          double *psi_2_real,
                                          double *psi_2_imag,
                                          double *error);