#define NVTX_USE true
#include "../../parallel_solver/crank_nicolson/cn_rect_psolver.cuh"
#include "../../../utils.h"
#include <iostream>
#include <cmath>

#include <string>
#include <fstream>
#include <thread>
#include <thread>
#include "nvToolsExt.h"

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
                                 float relaxation)
{
    __shared__ float tile_old_real[nTx][nTy];
    __shared__ float tile_old_imag[nTx][nTy];
    __shared__ float tile_new_real_trial[nTx][nTy];
    __shared__ float tile_new_imag_trial[nTx][nTy];
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
    tile_new_real_trial[thread_x][thread_y] = psi_new_real_trial[j * striding + i];
    tile_new_imag_trial[thread_x][thread_y] = psi_new_imag_trial[j * striding + i];
    tile_potential[thread_x][thread_y] = potential[j * striding + i];

    __syncthreads();

    // Parameters
    float sigma_x = tau / (4 * h_x * h_x);
    float sigma_y = tau / (4 * h_y * h_y);
    float a = 2 * sigma_x + 2 * sigma_y + 0.5 * tau * tile_potential[thread_x][thread_y];
    float b = a * a + 1;
    float amplitude_old = tile_old_real[thread_x][thread_y] * tile_old_real[thread_x][thread_y] + tile_old_imag[thread_x][thread_y] * tile_old_imag[thread_x][thread_y];
    float amplitude_new = tile_new_real_trial[thread_x][thread_y] * tile_new_real_trial[thread_x][thread_y] + tile_new_imag_trial[thread_x][thread_y] * tile_new_imag_trial[thread_x][thread_y];

    // Update tile position
    atomicAdd(&psi_new_real[j * striding + i],
              relaxation * ((1 - a * a) / b * tile_old_real[thread_x][thread_y] + 2 * a / b * tile_old_imag[thread_x][thread_y]) * (i >= 0) * (i < n_x) * (j >= 0) * (j < n_y));
    atomicAdd(&psi_new_real[j * striding + i],
              relaxation * (-a * g * tau / (2 * b) * tile_old_real[thread_x][thread_y] + g * tau / (2 * b) * tile_old_imag[thread_x][thread_y]) * amplitude_old * (i >= 0) * (i < n_x) * (j >= 0) * (j < n_y));
    atomicAdd(&psi_new_real[j * striding + i],
              relaxation * (-a * g * tau / (2 * b) * tile_new_real_trial[thread_x][thread_y] + g * tau / (2 * b) * tile_new_imag_trial[thread_x][thread_y]) * amplitude_new * (i >= 0) * (i < n_x) * (j >= 0) * (j < n_y));

    atomicAdd(&psi_new_imag[j * striding + i],
              relaxation * (-2 * a / b * tile_old_real[thread_x][thread_y] + (1 - a * a) / b * tile_old_imag[thread_x][thread_y]) * (i >= 0) * (i < n_x) * (j >= 0) * (j < n_y));
    atomicAdd(&psi_new_imag[j * striding + i],
              relaxation * (-g * tau / (2 * b) * tile_old_real[thread_x][thread_y] - a * g * tau / (2 * b) * tile_old_imag[thread_x][thread_y]) * amplitude_old * (i >= 0) * (i < n_x) * (j >= 0) * (j < n_y));
    atomicAdd(&psi_new_imag[j * striding + i],
              relaxation * (-g * tau / (2 * b) * tile_new_real_trial[thread_x][thread_y] - a * g * tau / (2 * b) * tile_new_imag_trial[thread_x][thread_y]) * amplitude_new * (i >= 0) * (i < n_x) * (j >= 0) * (j < n_y));

    // Update left
    atomicAdd(&psi_new_real[j * striding + i - 1],
              relaxation * (a * sigma_x / b * tile_old_real[thread_x][thread_y] - sigma_x / b * tile_old_imag[thread_x][thread_y]) *
                  (i > 0) * (i < n_x) * (j >= 0) * (j < n_y));
    atomicAdd(&psi_new_real[j * striding + i - 1],
              relaxation * (a * sigma_x / b * tile_new_real_trial[thread_x][thread_y] - sigma_x / b * tile_new_imag_trial[thread_x][thread_y]) *
                  (i > 0) * (i < n_x) * (j >= 0) * (j < n_y));

    atomicAdd(&psi_new_imag[j * striding + i - 1],
              relaxation * (sigma_x / b * tile_old_real[thread_x][thread_y] + a * sigma_x / b * tile_old_imag[thread_x][thread_y]) *
                  (i > 0) * (i < n_x) * (j >= 0) * (j < n_y));
    atomicAdd(&psi_new_imag[j * striding + i - 1],
              relaxation * (sigma_x / b * tile_new_real_trial[thread_x][thread_y] + a * sigma_x / b * tile_new_imag_trial[thread_x][thread_y]) *
                  (i > 0) * (i < n_x) * (j >= 0) * (j < n_y));

    // Update right
    atomicAdd(&psi_new_real[j * striding + i + 1],
              relaxation * (a * sigma_x / b * tile_old_real[thread_x][thread_y] - sigma_x / b * tile_old_imag[thread_x][thread_y]) *
                  (i >= 0) * (i < (n_x - 1)) * (j >= 0) * (j < n_y));
    atomicAdd(&psi_new_real[j * striding + i + 1],
              relaxation * (a * sigma_x / b * tile_new_real_trial[thread_x][thread_y] - sigma_x / b * tile_new_imag_trial[thread_x][thread_y]) *
                  (i >= 0) * (i < (n_x - 1)) * (j >= 0) * (j < n_y));

    atomicAdd(&psi_new_imag[j * striding + i + 1],
              relaxation * (sigma_x / b * tile_old_real[thread_x][thread_y] + a * sigma_x / b * tile_old_imag[thread_x][thread_y]) *
                  (i >= 0) * (i < (n_x - 1)) * (j >= 0) * (j < n_y));
    atomicAdd(&psi_new_imag[j * striding + i + 1],
              relaxation * (sigma_x / b * tile_new_real_trial[thread_x][thread_y] + a * sigma_x / b * tile_new_imag_trial[thread_x][thread_y]) *
                  (i >= 0) * (i < (n_x - 1)) * (j >= 0) * (j < n_y));

    // Update down
    atomicAdd(&psi_new_real[(j - 1) * striding + i],
              relaxation * (a * sigma_y / b * tile_old_real[thread_x][thread_y] - sigma_y / b * tile_old_imag[thread_x][thread_y]) *
                  (i >= 0) * (i < n_x) * (j > 0) * (j < n_y));
    atomicAdd(&psi_new_real[(j - 1) * striding + i],
              relaxation * (a * sigma_y / b * tile_new_real_trial[thread_x][thread_y] - sigma_y / b * tile_new_imag_trial[thread_x][thread_y]) *
                  (i >= 0) * (i < n_x) * (j > 0) * (j < n_y));

    atomicAdd(&psi_new_imag[(j - 1) * striding + i],
              relaxation * (sigma_y / b * tile_old_real[thread_x][thread_y] + a * sigma_y / b * tile_old_imag[thread_x][thread_y]) *
                  (i >= 0) * (i < n_x) * (j > 0) * (j < n_y));
    atomicAdd(&psi_new_imag[(j - 1) * striding + i],
              relaxation * (sigma_y / b * tile_new_real_trial[thread_x][thread_y] + a * sigma_y / b * tile_new_imag_trial[thread_x][thread_y]) *
                  (i >= 0) * (i < n_x) * (j > 0) * (j < n_y));

    // Update up
    atomicAdd(&psi_new_real[(j + 1) * striding + i],
              relaxation * (a * sigma_y / b * tile_old_real[thread_x][thread_y] - sigma_y / b * tile_old_imag[thread_x][thread_y]) *
                  (i >= 0) * (i < n_x) * (j >= 0) * (j < (n_y - 1)));
    atomicAdd(&psi_new_real[(j + 1) * striding + i],
              relaxation * (a * sigma_y / b * tile_new_real_trial[thread_x][thread_y] - sigma_y / b * tile_new_imag_trial[thread_x][thread_y]) *
                  (i >= 0) * (i < n_x) * (j >= 0) * (j < (n_y - 1)));

    atomicAdd(&psi_new_imag[(j + 1) * striding + i],
              relaxation * (sigma_y / b * tile_old_real[thread_x][thread_y] + a * sigma_y / b * tile_old_imag[thread_x][thread_y]) *
                  (i >= 0) * (i < n_x) * (j >= 0) * (j < (n_y - 1)));
    atomicAdd(&psi_new_imag[(j + 1) * striding + i],
              relaxation * (sigma_y / b * tile_new_real_trial[thread_x][thread_y] + a * sigma_y / b * tile_new_imag_trial[thread_x][thread_y]) *
                  (i >= 0) * (i < n_x) * (j >= 0) * (j < (n_y - 1)));
}
// Only works with single block
__global__ void calculate_local_error(float *psi_1_real,
                                      float *psi_1_imag,
                                      float *psi_2_real,
                                      float *psi_2_imag,
                                      float *error_array,
                                      int n_x,
                                      int n_y)
{

    __shared__ float tile_psi_1_real[nTx][nTy];
    __shared__ float tile_psi_1_imag[nTx][nTy];
    __shared__ float tile_psi_2_real[nTx][nTy];
    __shared__ float tile_psi_2_imag[nTx][nTy];

    int block_x = blockIdx.x * blockDim.x;
    int block_y = blockIdx.y * blockDim.y;
    int thread_x = threadIdx.x;
    int thread_y = threadIdx.y;
    int i = (block_x + thread_x);
    int j = (block_y + thread_y);
    int striding = gridDim.x * blockDim.x;

    tile_psi_1_real[thread_x][thread_y] = psi_1_real[j * striding + i];
    tile_psi_1_imag[thread_x][thread_y] = psi_1_imag[j * striding + i];
    tile_psi_2_real[thread_x][thread_y] = psi_2_real[j * striding + i];
    tile_psi_2_imag[thread_x][thread_y] = psi_2_imag[j * striding + i];

    __syncthreads();

    error_array[j * striding + i] = ((tile_psi_1_real[thread_x][thread_y] - tile_psi_2_real[thread_x][thread_y]) * (tile_psi_1_real[thread_x][thread_y] - tile_psi_2_real[thread_x][thread_y])) * (i >= 0) * (i < n_x) * (j >= 0) * (j < n_y);
    error_array[j * striding + i] += ((tile_psi_1_imag[thread_x][thread_y] - tile_psi_2_imag[thread_x][thread_y]) * (tile_psi_1_imag[thread_x][thread_y] - tile_psi_2_imag[thread_x][thread_y])) * (i >= 0) * (i < n_x) * (j >= 0) * (j < n_y);
}

__global__ void reduction_error(float *error_array, float *error, int array_size)
{

    int idx = threadIdx.x;
    float sum = 0;
    for (int i = idx; i < array_size; i += nTx * nTy)
    {
        sum += error_array[i];
    }

    __shared__ float r[nTx * nTy];
    r[idx] = sum;
    __syncthreads();
    for (int size = nTx * nTy / 2; size > 0; size /= 2)
    { // uniform
        if (idx < size)
            r[idx] += r[idx + size];
        __syncthreads();
    }
    if (idx == 0)
        *error = r[0];
}

__global__ void calculate_probability(float *psi_real, float *psi_imag, float *probability, int n_x, int n_y)
{

    __shared__ float tile_psi_real[nTx][nTy];
    __shared__ float tile_psi_imag[nTx][nTy];

    int block_x = blockIdx.x * blockDim.x;
    int block_y = blockIdx.y * blockDim.y;
    int thread_x = threadIdx.x;
    int thread_y = threadIdx.y;
    int i = (block_x + thread_x);
    int j = (block_y + thread_y);
    int striding = gridDim.x * blockDim.x;

    tile_psi_real[thread_x][thread_y] = psi_real[j * striding + i];
    tile_psi_imag[thread_x][thread_y] = psi_imag[j * striding + i];
    __syncthreads();

    probability[j * striding + i] = ((tile_psi_real[thread_x][thread_y] * tile_psi_real[thread_x][thread_y]) + (tile_psi_imag[thread_x][thread_y] * tile_psi_imag[thread_x][thread_y])) * (i >= 0) * (i < n_x) * (j >= 0) * (j < n_y);
}

// Only works with single block
__global__ void calculate_normalize_factor(float *probability, float *normalize_factor, int array_size, float unit_area)
{

    int idx = threadIdx.x;
    float sum = 0;
    for (auto i = idx; i < array_size; i += nTx * nTy)
    {
        sum += probability[i];
    }

    __shared__ float r[nTx * nTy];
    r[idx] = sum;
    __syncthreads();

    for (auto size = nTx * nTy / 2; size > 0; size /= 2)
    {
        if (idx < size)
            r[idx] += r[idx + size];
        __syncthreads();
    }
    if (idx == 0)
    {
        *normalize_factor = sqrt(r[0] * unit_area);
    }
}

__global__ void normalize(float *psi_real, float *psi_imag, float *normalize_factor)
{

    int block_x = blockIdx.x * blockDim.x;
    int block_y = blockIdx.y * blockDim.y;
    int thread_x = threadIdx.x;
    int thread_y = threadIdx.y;
    int i = (block_x + thread_x);
    int j = (block_y + thread_y);
    int striding = gridDim.x * blockDim.x;

    psi_real[j * striding + i] /= *normalize_factor;
    psi_imag[j * striding + i] /= *normalize_factor;
}

__global__ void scale_prev_solution(float *psi_real, float *psi_imag, float scale)
{

    int block_x = blockIdx.x * blockDim.x;
    int block_y = blockIdx.y * blockDim.y;
    int thread_x = threadIdx.x;
    int thread_y = threadIdx.y;
    int i = (block_x + thread_x);
    int j = (block_y + thread_y);
    int striding = gridDim.x * blockDim.x;

    psi_real[j * striding + i] *= scale;
    psi_imag[j * striding + i] *= scale;
}

CNRectPSolver::CNRectPSolver(
    // std::function<float(float, float)> potential,
    float g,
    RectangularDomain *domain_,
    int device_number)
    : BaseSolver(g)
{
    this->domain = domain_;
    this->string_info = std::string{"Crank_Nicolson_parallel_"};
    cudaSetDevice(device_number);
};

void fileout_debug(float *array, int n_x, int n_y, std::string filename)
{
    std::ofstream fileout(filename.data());
    for (auto i = 0; i < n_y; ++i)
    {
        for (auto j = 0; j < n_x - 1; ++j)
        {
            fileout << array[n_x * i + j] << ", ";
        }
        fileout << array[n_x * i + n_x - 1] << std::endl;
    }
}

// void CNRectPSolver::solve(float tolerance, int max_iter)
void CNRectPSolver::solve(float tolerance, int max_iter, std::string dir_name, bool print_info, bool save_data)
{
    cudaDeviceSynchronize();
    if (NVTX_USE)
    {
        nvtxRangePushA("solver started");
    }
    int n_x = this->domain->get_num_grid_1();
    int n_y = this->domain->get_num_grid_2();
    float dt = this->domain->get_dt();

    float *d_error;
    float *d_normalize_factor;
    float relaxation_parameter = 1.;

    dim3 TPB(nTx, nTy);
    dim3 nBlocks(n_x / nTx + (n_x % nTx != 0), n_y / nTy + (n_y % nTy != 0));

    float *d_psi_old_real, *d_psi_old_imag, *d_psi_new_real, *d_psi_new_real_trial, *d_psi_new_imag_trial, *d_psi_new_imag, *d_potential;
    float *d_probability_array, *d_error_array;

    if (NVTX_USE)
    {
        nvtxRangePushA("cuda malloc");
    }
    cudaMalloc((float **)&d_psi_old_real, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    cudaMalloc((float **)&d_psi_old_imag, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    cudaMalloc((float **)&d_psi_new_real_trial, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    cudaMalloc((float **)&d_psi_new_imag_trial, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    cudaMalloc((float **)&d_psi_new_real, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    cudaMalloc((float **)&d_psi_new_imag, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    cudaMalloc((float **)&d_potential, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    cudaMalloc((float **)&d_probability_array, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    cudaMalloc((float **)&d_error_array, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    cudaMalloc((float **)&d_error, sizeof(float));
    cudaMalloc((float **)&d_normalize_factor, sizeof(float));
    if (NVTX_USE)
    {
        nvtxRangePop();
    }

    cudaStream_t stream_psi_new_real, stream_psi_new_imag, stream_potential;
    cudaStream_t stream_device_to_device_1, stream_device_to_device_2, stream_device_to_device_3, stream_device_to_device_4;

    cudaStreamCreate(&stream_psi_new_real);
    cudaStreamCreate(&stream_psi_new_imag);
    cudaStreamCreate(&stream_potential);

    cudaStreamCreate(&stream_device_to_device_1);
    cudaStreamCreate(&stream_device_to_device_2);
    cudaStreamCreate(&stream_device_to_device_3);
    cudaStreamCreate(&stream_device_to_device_4);

    float *h_psi_new_real, *h_psi_new_imag, *h_potential;
    cudaHostAlloc((void **)&h_psi_new_real, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaHostAllocDefault);
    cudaHostAlloc((void **)&h_psi_new_imag, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaHostAllocDefault);
    cudaHostAlloc((void **)&h_potential, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaHostAllocDefault);

    if (NVTX_USE)
    {
        nvtxRangePushA("initialize: host array");
    }
    std::complex<float> wave_func;
    float potential_value;
    for (int i = 0; i < n_x; ++i)
    {
        for (int j = 0; j < n_y; ++j)
        {
            wave_func = this->domain->at(i, j, 0)->value;
            potential_value = this->domain->potential_grid->at(i, j)->value.real();
            h_psi_new_real[j * TPB.x * nBlocks.x + i] = wave_func.real();
            h_psi_new_imag[j * TPB.x * nBlocks.x + i] = wave_func.imag();
            h_potential[j * TPB.x * nBlocks.x + i] = potential_value;
        }
    }

    if (NVTX_USE)
    {
        nvtxRangePop();
    }

    if (NVTX_USE)
    {
        nvtxRangePushA("memcpyAsync");
    }
    cudaMemcpyAsync(d_psi_new_real, h_psi_new_real, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyHostToDevice, stream_psi_new_real);
    cudaMemcpyAsync(d_psi_new_imag, h_psi_new_imag, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyHostToDevice, stream_psi_new_imag);
    cudaMemcpyAsync(d_potential, h_potential, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyHostToDevice, stream_potential);
    cudaMemcpyAsync(d_psi_new_real_trial, d_psi_new_real, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyDeviceToDevice, stream_device_to_device_1);
    cudaMemcpyAsync(d_psi_new_imag_trial, d_psi_new_imag, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyDeviceToDevice, stream_device_to_device_2);
    cudaMemcpyAsync(d_psi_old_real, d_psi_new_real, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyDeviceToDevice, stream_device_to_device_3);
    cudaMemcpyAsync(d_psi_old_imag, d_psi_new_imag, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyDeviceToDevice, stream_device_to_device_4);

    if (NVTX_USE)
    {
        nvtxRangePop();
    }

    if (save_data)
    {
        this->domain->generate_directory_name(this->string_info + dir_name, print_info);
        // Save initial condition
        this->domain->generate_single_txt_file(std::string("Solution_") + std::to_string(0));
    }
    else
    {
        this->domain->update_time(true);
    }
    float *buffer_real = (float *)malloc(sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    float *buffer_imag = (float *)malloc(sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);

    ///////////////////////////////////////  Calculation Core started ////////////////////////////////////////////////////////////
    std::vector<std::thread> threads;
    threads.reserve(this->domain->get_num_times() - 1);
    for (auto k = 0; k < this->domain->get_num_times() - 1; ++k)
    {
        {
            this->solve_single_time(k, d_psi_old_real,
                                    d_psi_old_imag,
                                    d_psi_new_real_trial,
                                    d_psi_new_imag_trial,
                                    d_psi_new_real,
                                    d_psi_new_imag,
                                    d_potential,
                                    d_probability_array,
                                    d_normalize_factor,
                                    d_error_array,
                                    d_error,
                                    max_iter,
                                    tolerance,
                                    relaxation_parameter,
                                    nBlocks,
                                    TPB,
                                    stream_device_to_device_1,
                                    stream_device_to_device_2,
                                    buffer_real,
                                    buffer_imag,
                                    save_data);

            auto export_thread = std::thread(&CNRectPSolver::export_single_time,
                                             this,
                                             k,
                                             buffer_real,
                                             buffer_imag,
                                             nBlocks,
                                             TPB,
                                             save_data);
            threads.push_back(std::move(export_thread));
        }
    }

    for (std::thread &th : threads)
    {
        if (th.joinable())
            th.join();
    }

    ///////////////////////////////////////  Calculation Core finished ////////////////////////////////////////////////////////////
    cudaFreeHost(h_psi_new_real);
    cudaFreeHost(h_psi_new_imag);
    cudaFreeHost(h_potential);

    if (NVTX_USE)
    {
        nvtxRangePushA("cuda free");
    }

    cudaFree(d_psi_old_real);
    cudaFree(d_psi_old_imag);
    cudaFree(d_psi_new_real_trial);
    cudaFree(d_psi_new_imag_trial);
    cudaFree(d_psi_new_real);
    cudaFree(d_psi_new_imag);
    cudaFree(d_potential);
    cudaFree(d_probability_array);
    cudaFree(d_error_array);
    cudaFree(d_error);
    cudaFree(d_normalize_factor);

    if (NVTX_USE)
    {
        nvtxRangePop();
    }

    if (print_info)
    {
        this->domain->print_directory_info();
    }
    if (NVTX_USE)
    {
        nvtxRangePop();
    }
}

void CNRectPSolver::solve_single_time(int k,
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
                                      bool save_data)
{
    int n_x = this->domain->get_num_grid_1();
    int n_y = this->domain->get_num_grid_2();

    float h_x = this->domain->get_infinitesimal_distance1();
    float h_y = this->domain->get_infinitesimal_distance2();

    if (k > 0)
    {
        if (NVTX_USE)
        {
            nvtxRangePushA("memcpy new to old");
        }
        cudaMemcpyAsync(d_psi_old_real, d_psi_new_real, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyDeviceToDevice, stream_device_to_device_1);
        cudaMemcpyAsync(d_psi_old_imag, d_psi_new_imag, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyDeviceToDevice, stream_device_to_device_2);
        if (NVTX_USE)
        {
            nvtxRangePop();
        }
    }
    if (NVTX_USE)
    {
        nvtxRangePushA((std::string("solve time ") + std::to_string(k)).c_str());
    }

    float error = 1.;
    for (auto iter = 0; iter < max_iter; ++iter)
    {
        if (error < tolerance)
        {
            break;
        }
        if (NVTX_USE)
        {
            nvtxRangePushA("memcpy new to trial");
        }
        cudaMemcpyAsync(d_psi_new_real_trial, d_psi_new_real, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyDeviceToDevice, stream_device_to_device_1);
        cudaMemcpyAsync(d_psi_new_imag_trial, d_psi_new_imag, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyDeviceToDevice, stream_device_to_device_2);
        if (NVTX_USE)
        {
            nvtxRangePop();
        }

        if (NVTX_USE)
        {
            nvtxRangePushA("scale_prev_solution");
        }
        scale_prev_solution<<<nBlocks, TPB>>>(d_psi_new_real, d_psi_new_imag, 1 - relaxation_parameter);
        if (NVTX_USE)
        {
            nvtxRangePop();
        }

        cudaDeviceSynchronize();
        if (NVTX_USE)
        {
            nvtxRangePushA("cn_rect_cusolver");
        }
        cn_rect_cusolver<<<nBlocks, TPB>>>(
            d_psi_old_real,
            d_psi_old_imag,
            d_psi_new_real_trial,
            d_psi_new_imag_trial,
            d_psi_new_real,
            d_psi_new_imag,
            d_potential,
            n_x, n_y,
            this->g,
            this->domain->get_infinitesimal_distance1(),
            this->domain->get_infinitesimal_distance2(),
            this->domain->get_dt(),
            relaxation_parameter);
        if (NVTX_USE)
        {
            nvtxRangePop();
        }
        if (NVTX_USE)
        {
            nvtxRangePushA("initialize prob and normalize factor");
        }
        cudaMemset(d_probability_array, 0, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
        cudaMemset(d_normalize_factor, 0, sizeof(float));
        if (NVTX_USE)
        {
            nvtxRangePop();
        }

        if (NVTX_USE)
        {
            nvtxRangePushA("calculate_probability");
        }
        calculate_probability<<<nBlocks, TPB>>>(d_psi_new_real, d_psi_new_imag, d_probability_array, n_x, n_y);
        if (NVTX_USE)
        {
            nvtxRangePop();
        }

        if (NVTX_USE)
        {
            nvtxRangePushA("calculate_normalize_factor");
        }
        calculate_normalize_factor<<<1, TPB.x * TPB.y>>>(d_probability_array, d_normalize_factor, TPB.x * nBlocks.x * TPB.y * nBlocks.y, h_x * h_y);
        if (NVTX_USE)
        {
            nvtxRangePop();
        }

        if (NVTX_USE)
        {
            nvtxRangePushA("normalize");
        }
        normalize<<<nBlocks, TPB>>>(d_psi_new_real, d_psi_new_imag, d_normalize_factor);
        if (NVTX_USE)
        {
            nvtxRangePop();
        }

        cudaDeviceSynchronize();

        if (NVTX_USE)
        {
            nvtxRangePushA("initialize local and global error");
        }
        cudaMemset(d_error_array, 0, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
        cudaMemset(d_error, 0, sizeof(float));
        if (NVTX_USE)
        {
            nvtxRangePop();
        }

        if (NVTX_USE)
        {
            nvtxRangePushA("calculate_local_error");
        }
        calculate_local_error<<<nBlocks, TPB>>>(d_psi_new_real, d_psi_new_imag, d_psi_new_real_trial, d_psi_new_imag_trial, d_error_array, n_x, n_y);
        if (NVTX_USE)
        {
            nvtxRangePop();
        }

        if (NVTX_USE)
        {
            nvtxRangePushA("reduction_error");
        }
        reduction_error<<<1, TPB.x * TPB.y>>>(d_error_array, d_error, TPB.x * nBlocks.x * TPB.y * nBlocks.y);
        if (NVTX_USE)
        {
            nvtxRangePop();
        }

        cudaMemcpy(&error, d_error, sizeof(float), cudaMemcpyDeviceToHost);
    }
    cudaMemcpyAsync(buffer_real, d_psi_new_real, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyDeviceToHost, stream_device_to_device_1);
    cudaMemcpyAsync(buffer_imag, d_psi_new_imag, sizeof(float) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyDeviceToHost, stream_device_to_device_2);

    if (NVTX_USE)
    {
        nvtxRangePop();
    }
}

void CNRectPSolver::export_single_time(int k,
                                       float *buffer_real,
                                       float *buffer_imag,
                                       dim3 nBlocks,
                                       dim3 TPB,
                                       bool save_data)
{
    int n_x = this->domain->get_num_grid_1();
    int n_y = this->domain->get_num_grid_2();

    if (NVTX_USE)
    {
        nvtxRangePushA("save_final_data");
    }
    for (int i = 0; i < n_x; ++i)
    {

        for (int j = 0; j < n_y; ++j)
        {
            this->domain->assign_wave_function(i, j, k + 1,
                                               std::complex<float>{buffer_real[j * TPB.x * nBlocks.x + i],
                                                                   buffer_imag[j * TPB.x * nBlocks.x + i]});
        }
    }
    if (NVTX_USE)
    {
        nvtxRangePop();
    }
    if (NVTX_USE)
    {
        nvtxRangePushA("update time");
    }
    if (save_data)
    {
        this->domain->generate_single_txt_file(std::string("Solution_") + std::to_string(k + 1), true);
    }

    this->domain->update_time(true);

    if (NVTX_USE)
    {
        nvtxRangePop();
    }
}