#pragma once
#include <complex>
#include <vector>
#include <functional>
#include <iostream>
#include <cmath>

#include "../../../domain/rect_domain.h"
#include "../../../initial_condition/initial_condition.h"
#include "../../base_solver.h"
#include "../../serial_solver/forward_euler/fe_rect_solver.h"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#define nTx 16
#define nTy 16

class CNRectPSolver : BaseSolver
{
public:
    CNRectPSolver() = default;
    CNRectPSolver(
        // std::function<float(float, float)> potential,
        float g,
        RectangularDomain *domain);
    void generateRectangularDomain();
    void solve(float tolerance, int max_iter);

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

__global__ void cn_rect_cusolver(float *psi_old_real,
                                 float *psi_old_imag,
                                 float *psi_new_real_trial,
                                 float *psi_new_imag_trial,
                                 float *psi_new_real,
                                 float *psi_new_imag,
                                 float *potential,
                                 int n_x,
                                 int n_y,
                                 float g,
                                 float h_x,
                                 float h_y,
                                 float tau,
                                 float relaxation);

__global__ void calculate_error_on_device(float *psi_1_real,
                                          float *psi_1_imag,
                                          float *psi_2_real,
                                          float *psi_2_imag,
                                          float *error, int n_x, int n_y);

__global__ void calculate_local_error(float *psi_1_real,
                                      float *psi_1_imag,
                                      float *psi_2_real,
                                      float *psi_2_imag,
                                      float *error_array,
                                      int n_x,
                                      int n_y);

__global__ void reduction_error(float *error_array, float *error, int array_size);

__global__ void calculate_probability(float *psi_real, float *psi_imag, float *probability, int n_x, int n_y);

__global__ void calculate_normalize_factor(float *probability, float *normalize_factor, int array_size, float unit_area);

__global__ void normalize(float *psi_real, float *psi_imag, float *normalize_factor);

__global__ void scale_prev_solution(float *psi_real, float *psi_imag, float scale);

void fileout_debug(float *array, int n_x, int n_y, std::string filename);