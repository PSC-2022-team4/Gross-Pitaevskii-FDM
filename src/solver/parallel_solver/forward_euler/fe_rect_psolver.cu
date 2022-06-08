#include "fe_rect_psolver.cuh"
#include <iostream>
#include <string>
#include <stdio.h>
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
                                 float tau)
{
    __shared__ float tile_old_real[nTx][nTy];
    __shared__ float tile_old_imag[nTx][nTy];
    __shared__ float tile_potential[nTx][nTy];

    int block_x = blockIdx.x * blockDim.x;
    int block_y = blockIdx.y * blockDim.y;
    int thread_x = threadIdx.x;
    int thread_y = threadIdx.y;
    int i = (block_x + thread_x);
    int j = (block_y + thread_y);
    int striding = gridDim.x * blockDim.x;

    tile_old_real[thread_x][thread_y] = psi_old_real[j * striding + i];
    tile_old_imag[thread_x][thread_y] = psi_old_imag[j * striding + i];
    tile_potential[thread_x][thread_y] = potential[j * striding + i];

    __syncthreads();

    //  Update tile position
    psi_new_real[j * striding + i] += ((1 -
                                        tile_potential[thread_x][thread_y] -
                                        g * (tile_old_real[thread_x][thread_y] * tile_old_real[thread_x][thread_y] + tile_old_imag[thread_x][thread_y] * tile_old_imag[thread_x][thread_y])) *
                                           tile_old_real[thread_x][thread_y] +
                                       (2 * tau / (h_x * h_x) + 2 * tau / (h_y * h_y)) * tile_old_imag[thread_x][thread_y]) *
                                      (i >= 0) * (i < n_x) * (j >= 0) * (j < n_y);
    psi_new_imag[j * striding + i] += ((1 -
                                        tile_potential[thread_x][thread_y] -
                                        g * (tile_old_real[thread_x][thread_y] * tile_old_real[thread_x][thread_y] + tile_old_imag[thread_x][thread_y] * tile_old_imag[thread_x][thread_y])) *
                                           tile_old_imag[thread_x][thread_y] -
                                       (2 * tau / (h_x * h_x) + 2 * tau / (h_y * h_y)) * tile_old_real[thread_x][thread_y]) *
                                      (i >= 0) * (i < n_x) * (j >= 0) * (j < n_y);

    // Update left position
    psi_new_real[j * striding + i - 1] += (-(tau / (h_x * h_x)) * tile_old_imag[thread_x][thread_y]) * (i > 0) * (i < n_x) * (j >= 0) * (j < n_y);
    psi_new_imag[j * striding + i - 1] += ((tau / (h_x * h_x)) * tile_old_real[thread_x][thread_y]) * (i > 0) * (i < n_x) * (j >= 0) * (j < n_y);

    // Update right position
    psi_new_real[j * striding + i + 1] += (-(tau / (h_x * h_x)) * tile_old_imag[thread_x][thread_y]) * (i >= 0) * (i < (n_x - 1)) * (j >= 0) * (j < n_y);
    psi_new_imag[j * striding + i + 1] += ((tau / (h_x * h_x)) * tile_old_real[thread_x][thread_y]) * (i >= 0) * (i < (n_x - 1)) * (j >= 0) * (j < n_y);

    // Update down position
    psi_new_real[(j - 1) * striding + i] += (-(tau / (h_y * h_y)) * tile_old_imag[thread_x][thread_y]) * (i >= 0) * (i < n_x) * (j > 0) * (j < n_y);
    psi_new_imag[(j - 1) * striding + i] += ((tau / (h_y * h_y)) * tile_old_real[thread_x][thread_y]) * (i >= 0) * (i < n_x) * (j > 0) * (j < n_y);

    // Update up position
    psi_new_real[(j + 1) * striding + i] += (-(tau / (h_y * h_y)) * tile_old_imag[thread_x][thread_y]) * (i >= 0) * (i < n_x) * (j >= 0) * (j < (n_y - 1));
    psi_new_imag[(j + 1) * striding + i] += ((tau / (h_y * h_y)) * tile_old_real[thread_x][thread_y]) * (i >= 0) * (i < n_x) * (j >= 0) * (j < (n_y - 1));
}

FERectPSolver::FERectPSolver(
    // std::function<float(float, float)> potential_,
    float g_,
    RectangularDomain *domain_,
    int device_number
    )
    : BaseSolver(g_)
{
    this->domain = domain_;
    // this->generate_potential_grid();
};
// void FERectPSolver::generate_potential_grid()
// {
//     int num_grid_1 = this->domain->get_num_grid_1();
//     int num_grid_2 = this->domain->get_num_grid_2();
//     float x_start = this->domain->at(0, 0, 0)->x;
//     float y_start = this->domain->at(0, 0, 0)->y;
//     float x_end = this->domain->at(num_grid_1 - 1, num_grid_2 - 1, 0)->x;
//     float y_end = this->domain->at(num_grid_1 - 1, num_grid_2 - 1, 0)->y;
//     this->potential_grid = RectangularSpatialGrid(num_grid_1, num_grid_2, x_start, x_end, y_start, y_end);
//     for (auto i = 0; i < num_grid_1; ++i)
//     {
//         for (auto j = 0; j < num_grid_2; ++j)
//         {
//             auto point = potential_grid.at(i, j);
//             point->value = {this->potential_func(point->x, point->y), 0};
//         }
//     }
// };
float FERectPSolver::get_potential_value(int i, int j)
{
    return this->domain->potential_grid.at(i, j)->value.real();
}


void FERectPSolver::solve_single_time(int k)
{
    int n_x = this->domain->get_num_grid_1();
    int n_y = this->domain->get_num_grid_2();
    float dt = this->domain->get_dt();

    dim3 TPB(nTx, nTy);
    dim3 nBlocks(n_x / nTx + (n_x % nTx != 0), n_y / nTy + (n_y % nTy != 0));
    float *h_psi_old_real = (float *)malloc(sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    float *h_psi_old_imag = (float *)malloc(sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    float *h_psi_new_real = (float *)malloc(sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    float *h_psi_new_imag = (float *)malloc(sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    float *h_potential = (float *)malloc(sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    float *d_psi_old_real, *d_psi_old_imag, *d_psi_new_real, *d_psi_new_imag, *d_potential;
    cudaMalloc((float **)&d_psi_old_real, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    cudaMalloc((float **)&d_psi_old_imag, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    cudaMalloc((float **)&d_psi_new_real, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    cudaMalloc((float **)&d_psi_new_imag, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    cudaMalloc((float **)&d_potential, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);

    std::complex<float> wave_func;
    float potential_value;
    for (int i = 0; i < n_x; ++i)
    {
        for (int j = 0; j < n_y; ++j)
        {
            wave_func = this->domain->at(i, j, k)->value;
            potential_value = this->domain->potential_grid.at(i, j)->value.real();
            h_psi_old_real[j * TPB.x * nBlocks.x + i] = wave_func.real();
            h_psi_old_imag[j * TPB.x * nBlocks.x + i] = wave_func.imag();
            h_psi_new_real[j * TPB.x * nBlocks.x + i] = 0.;
            h_psi_new_imag[j * TPB.x * nBlocks.x + i] = 0.;
            h_potential[j * TPB.x * nBlocks.x + i] = potential_value;
        }
    }

    cudaMemcpy(d_psi_new_real, h_psi_new_real, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyHostToDevice);
    cudaMemcpy(d_psi_new_imag, h_psi_new_imag, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyHostToDevice);
    cudaMemcpy(d_psi_old_real, h_psi_old_real, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyHostToDevice);
    cudaMemcpy(d_psi_old_imag, h_psi_old_imag, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyHostToDevice);
    cudaMemcpy(d_potential, h_potential, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyHostToDevice);
    cudaDeviceSynchronize();
    fe_rect_cusolver<<<nBlocks, TPB>>>(
        d_psi_old_real,
        d_psi_old_imag,
        d_psi_new_real,
        d_psi_new_imag,
        d_potential,
        n_x, n_y,
        this->g,
        this->domain->get_infinitesimal_distance1(),
        this->domain->get_infinitesimal_distance2(),
        this->domain->get_dt());
    cudaDeviceSynchronize();
    cudaMemcpy(h_psi_new_real, d_psi_new_real, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_psi_new_imag, d_psi_new_imag, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyDeviceToHost);

    for (int i = 0; i < n_x; ++i)
    {
        for (int j = 0; j < n_y; ++j)
        {
            this->domain->at(i, j, k + 1)->value = h_psi_new_real[j * TPB.x * nBlocks.x + i] + std::complex<float>{0, 1.} * h_psi_new_imag[j * TPB.x * nBlocks.x + i];
        }
    }
}

void FERectPSolver::solve(std::string dir_name)
{
    int time_length = this->domain->get_num_times();

    for (int k = 0; k < time_length - 1; ++k)
    {
        // std::cout << "Time step: " << k << std::endl;
        this->solve_single_time(k);
        this->domain->normalize(k + 1);
        this->domain->update_time();
    }
    // this->domain->generate_txt_file(std::string{"Forward_Euler_Result"} + dir_name);
}
