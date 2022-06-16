/**
 * @file cn_rect_psolver.cuh
 * @author Gyeonghun Kim, Minyoung Kim
 * @brief Header file for CUDA based parallel crank nicolson solver
 * @version 0.1
 * @date 2022-06-05
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#pragma once
#include <complex>
#include <vector>
#include <functional>
#include <iostream>
#include <cmath>

#include <string>
#include <fstream>
#include <thread>

#include "nvToolsExt.h"

#include "../../../utils.h"
#include "../../../domain/rect_domain.h"
#include "../../../initial_condition/initial_condition.h"
#include "../../base_solver.h"
#include "../../serial_solver/forward_euler/fe_rect_solver.h"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#define nTx 16
#define nTy 16

class CNRectPSolver : public BaseSolver
{
public:
    CNRectPSolver() = default;
    CNRectPSolver(
        // std::function<float(float, float)> potential,
        float g,
        RectangularDomain *domain,
        int device_number);
    // void generateRectangularDomain();
    // void solve(float tolerance, int max_iter);
    void solve(float tolerance, int max_iter, std::string dir_name = "", bool print_info = true, bool save_data = true);
    void solve_single_time(int time_index,
                           float *d_psi_old_real,
                           float *d_psi_old_imag,
                           float *d_psi_new_real_trial,
                           float *d_psi_new_imag_trial,
                           float *d_psi_new_real,
                           float *d_psi_new_imag,
                           float *d_potential,
                           float *d_probability_array,
                           float *d_normalize_factor,
                           float *d_error_array,
                           float *d_error,
                           int max_iter,
                           double tolerance,
                           double relaxation_parameter,
                           dim3 nBlocks,
                           dim3 TPB,
                           cudaStream_t stream_device_to_device_1,
                           cudaStream_t stream_device_to_device_2,
                           float *buffer_real,
                           float *buffer_imag,
                           bool save_data);

    void export_single_time(int k,
                            float *buffer_real,
                            float *buffer_imag,
                            dim3 nBlocks,
                            dim3 TPB,
                            bool save_data);
    // void solve_single_time(int k, float tolerance, int max_iter);

protected:
    RectangularDomain *domain;
    // RectangularSpatialGrid potential_grid;
    RectangularSpatialGrid *guess;
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