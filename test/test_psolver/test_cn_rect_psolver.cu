#include "test_cn_rect_psolver.cuh"
#include "../../src/solver/parallel_solver/crank_nicolson/cn_rect_psolver.cuh"
#include "../../src/potential/harmonic_potential.h"
#include "../../src/utils.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>

#include <functional>
#include <iostream>
#include <complex>
#include "gtest/gtest.h"

TEST(CNPSolverTest, InitializeSolveTest)
{
    bool all_passed = true;
    RectangularDomain *domain = (new RectangularDomain(32, 32, 0, 1, 3, -10, 10, -10, 10));
    auto initial_cond_function = [](float x, float y)
    { return std::complex<float>{1. * expf(-(x * x + y * y) / (1))}; };
    auto *initial_condition = new InitialCondition(initial_cond_function);
    initial_condition->assign_to_domain(domain);

    auto *potential = new HarmonicPotential(3, 5);
    potential->calcualte_potential_in_grid(domain);

    float g = -1;
    CNRectPSolver solver = CNRectPSolver(g, domain, 0);

    solver.solve(1e-11, 101);
    ASSERT_TRUE(all_passed);
}

TEST(CNPSolverTest, NormalizeTest)
{
    int n_x = 17;
    int n_y = 17;
    float h_x = 1;
    float h_y = 1;

    dim3 TPB(nTx, nTy);
    dim3 nBlocks(n_x / nTx + (n_x % nTx != 0), n_y / nTy + (n_y % nTy != 0));
    float *real_array = (float *)malloc(sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    float *imag_array = (float *)malloc(sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    float *prob_array = (float *)malloc(sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    float normalize_factor = 0;

    for (auto i = 0; i < n_x; ++i)
    {
        for (auto j = 0; j < n_y; ++j)
        {
            real_array[TPB.x * nBlocks.x * j + i] = 1.;
            imag_array[TPB.x * nBlocks.x * j + i] = 1.;
        }
    }
    float *d_real_array, *d_imag_array, *d_prob_array, *d_normalize_factor;

    cudaMalloc((float **)&d_real_array, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    cudaMalloc((float **)&d_imag_array, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    cudaMalloc((float **)&d_prob_array, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    cudaMalloc((float **)&d_normalize_factor, sizeof(float));
    cudaMemcpy(d_real_array, real_array, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyHostToDevice);
    cudaMemcpy(d_imag_array, imag_array, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyHostToDevice);
    cudaMemset(&d_prob_array, 0, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    cudaMemset(&d_normalize_factor, 0, sizeof(float));
    calculate_probability<<<nBlocks, TPB>>>(d_real_array, d_imag_array, d_prob_array, n_x, n_y);
    calculate_normalize_factor<<<1, TPB.x * TPB.y>>>(d_prob_array, d_normalize_factor, TPB.x * nBlocks.x * TPB.y * nBlocks.y, h_x * h_y);
    cudaMemcpy(&normalize_factor, d_normalize_factor, sizeof(float), cudaMemcpyDeviceToHost);

    ASSERT_FLOAT_EQ(normalize_factor, sqrt(2 * n_x * n_y));

    normalize<<<nBlocks, TPB>>>(d_real_array, d_imag_array, d_normalize_factor);
    cudaMemcpy(real_array, d_real_array, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyDeviceToHost);
    cudaMemcpy(imag_array, d_imag_array, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyDeviceToHost);

    ASSERT_FLOAT_EQ(real_array[0], 1 / sqrt(2 * n_x * n_y));
    ASSERT_FLOAT_EQ(imag_array[0], 1 / sqrt(2 * n_x * n_y));
    ASSERT_TRUE(abs(real_array[n_x + 1]) < 1e-6);

    cudaFree(d_real_array);
    cudaFree(d_imag_array);
    cudaFree(d_normalize_factor);
    cudaFree(d_prob_array);
}
TEST(CNPSolverTest, ErrorEstimation)
{
    int n_x = 17;
    int n_y = 17;

    dim3 TPB(nTx, nTy);
    dim3 nBlocks(n_x / nTx + (n_x % nTx != 0), n_y / nTy + (n_y % nTy != 0));
    float *real_array_1 = (float *)malloc(sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    float *imag_array_1 = (float *)malloc(sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    float *real_array_2 = (float *)malloc(sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    float *imag_array_2 = (float *)malloc(sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    float *error_array = (float *)malloc(sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    float error;
    for (auto i = 0; i < n_x; ++i)
    {
        for (auto j = 0; j < n_y; ++j)
        {
            real_array_1[TPB.x * nBlocks.x * j + i] = 1.;
            imag_array_1[TPB.x * nBlocks.x * j + i] = 1.;
            real_array_2[TPB.x * nBlocks.x * j + i] = 2.;
            imag_array_2[TPB.x * nBlocks.x * j + i] = 2.;
        }
    }
    float *d_real_array_1, *d_imag_array_1, *d_real_array_2, *d_imag_array_2, *d_error_array, *d_error;

    cudaMalloc((float **)&d_real_array_1, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    cudaMalloc((float **)&d_imag_array_1, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    cudaMalloc((float **)&d_real_array_2, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    cudaMalloc((float **)&d_imag_array_2, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    cudaMalloc((float **)&d_error_array, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    cudaMalloc((float **)&d_error, sizeof(float));

    cudaMemcpy(d_real_array_1, real_array_1, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyHostToDevice);
    cudaMemcpy(d_imag_array_1, imag_array_1, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyHostToDevice);
    cudaMemcpy(d_real_array_2, real_array_2, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyHostToDevice);
    cudaMemcpy(d_imag_array_2, imag_array_2, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyHostToDevice);
    cudaMemset(&d_error_array, 0, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    cudaMemset(&d_error, 0, sizeof(float));

    calculate_local_error<<<nBlocks, TPB>>>(d_real_array_1, d_imag_array_1, d_real_array_2, d_imag_array_2, d_error_array, n_x, n_y);
    cudaMemcpy(error_array, d_error_array, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyDeviceToHost);
    ASSERT_FLOAT_EQ(error_array[0], 2);
    ASSERT_FLOAT_EQ(error_array[n_x], 0);
    reduction_error<<<1, TPB.x * TPB.y>>>(d_error_array, d_error, TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    cudaMemcpy(&error, d_error, sizeof(float), cudaMemcpyDeviceToHost);
    ASSERT_FLOAT_EQ(error, 2 * n_x * n_y);

    cudaFree(d_real_array_1);
    cudaFree(d_imag_array_1);
    cudaFree(d_real_array_2);
    cudaFree(d_imag_array_2);
    cudaFree(d_error_array);
    cudaFree(d_error);
}