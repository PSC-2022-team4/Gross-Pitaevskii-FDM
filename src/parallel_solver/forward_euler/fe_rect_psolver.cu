#include "fe_rect_psolver.cuh"
#include <iostream>
#include <string>
#include <stdio.h>
__global__ void fe_rect_cusolver(double *psi_old_real,
                                 double *psi_old_imag,
                                 double *psi_new_real,
                                 double *psi_new_imag,
                                 double *potential,
                                 int n_x,
                                 int n_y,
                                 double g,
                                 double h_x,
                                 double h_y,
                                 double tau)
{
    __shared__ double tile_old_real[nTx][nTy];
    __shared__ double tile_old_imag[nTx][nTy];
    __shared__ double tile_potential[nTx][nTy];

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
    std::function<double(double, double)> potential_,
    double g_,
    RectangularDomain *domain_)
    : BaseSolver(potential_, g_)
{
    this->domain = domain_;
    this->generate_potential_grid();
};
void FERectPSolver::generate_potential_grid()
{
    int num_grid_1 = this->domain->get_num_grid_1();
    int num_grid_2 = this->domain->get_num_grid_2();
    double x_start = this->domain->at(0, 0, 0)->x;
    double y_start = this->domain->at(0, 0, 0)->y;
    double x_end = this->domain->at(num_grid_1 - 1, num_grid_2 - 1, 0)->x;
    double y_end = this->domain->at(num_grid_1 - 1, num_grid_2 - 1, 0)->y;
    this->potential_grid = RectangularSpatialGrid(num_grid_1, num_grid_2, x_start, x_end, y_start, y_end);
    for (auto i = 0; i < num_grid_1; ++i)
    {
        for (auto j = 0; j < num_grid_2; ++j)
        {
            auto point = potential_grid.at(i, j);
            point->wave_function = {this->potential_func(point->x, point->y), 0};
        }
    }
};
double FERectPSolver::get_potential_value(int i, int j)
{
    return this->potential_grid.at(i, j)->wave_function.real();
}
/**
 * @brief Time differential of phi
 *
 * @param i index for x
 * @param j index for y
 * @param k index for time(t)
 * @return std::complex<double> time differential at x, y, t
 */
std::complex<double> FERectPSolver::temporal_equation(int i, int j, int k)
{
    // Use five stencil method
    auto point_data = this->domain->at(i, j, k);

    // l,r,d,u denotes left, right, down, up value
    // Check boundary
    auto point_data_l = this->domain->at(i - 1, j, k);
    if (i <= 0)
        point_data_l = new GridPoint(0., 0., std::complex<double>{0, 0});
    auto point_data_d = this->domain->at(i, j - 1, k);
    if (j <= 0)
        point_data_d = new GridPoint(0., 0., std::complex<double>{0, 0});
    auto point_data_r = this->domain->at(i + 1, j, k);
    if (i >= (this->domain->get_num_grid_1()) - 1)
        point_data_r = new GridPoint(0., 0., std::complex<double>{0, 0});
    auto point_data_u = this->domain->at(i, j + 1, k);
    if (j >= (this->domain->get_num_grid_2()) - 1)
        point_data_u = new GridPoint(0., 0., std::complex<double>{0, 0});

    // potential at x, y
    double V_ij = this->get_potential_value(i, j);
    // this->potential_func(point_data->x, point_data->y);

    // g * |psi(x,y)|^2
    double additional_term = (this->g) * (std::abs(point_data->wave_function)) * (std::abs(point_data->wave_function));

    // Set infinitesimal value
    double dx = this->domain->get_infinitesimal_distance1();
    double dy = this->domain->get_infinitesimal_distance2();
    // df denote time differential of dt (d(psi)/dt)
    //  = (laplace - V-g|psi|^2) psi
    std::complex<double> df =
        +((point_data_r->wave_function) + (point_data_l->wave_function) - (point_data->wave_function) * std::complex<double>{2}) / (std::complex<double>{dx * dx}) + ((point_data_u->wave_function) + (point_data_d->wave_function) - (point_data->wave_function) * std::complex<double>{2}) / (std::complex<double>{dy * dy}) - (V_ij + additional_term) * (point_data->wave_function);
    df *= std::complex<double>{0, 1};
    return df;
};

void FERectPSolver::solve_single_time(int k)
{
    int n_x = this->domain->get_num_grid_1();
    int n_y = this->domain->get_num_grid_2();
    double dt = this->domain->get_dt();

    dim3 TPB(nTx, nTy);
    dim3 nBlocks(n_x / nTx + (n_x % nTx != 0), n_y / nTy + (n_y % nTy != 0));
    double *h_psi_old_real = (double *)malloc(sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    double *h_psi_old_imag = (double *)malloc(sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    double *h_psi_new_real = (double *)malloc(sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    double *h_psi_new_imag = (double *)malloc(sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    double *h_potential = (double *)malloc(sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    double *d_psi_old_real, *d_psi_old_imag, *d_psi_new_real, *d_psi_new_imag, *d_potential;
    cudaMalloc((double **)&d_psi_old_real, sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    cudaMalloc((double **)&d_psi_old_imag, sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    cudaMalloc((double **)&d_psi_new_real, sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    cudaMalloc((double **)&d_psi_new_imag, sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    cudaMalloc((double **)&d_potential, sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);

    std::complex<double> wave_func;
    double potential_value;
    for (int i = 0; i < n_x; ++i)
    {
        for (int j = 0; j < n_y; ++j)
        {
            wave_func = this->domain->at(i, j, k)->wave_function;
            potential_value = this->potential_grid.at(i, j)->wave_function.real();
            h_psi_old_real[j * TPB.x * nBlocks.x + i] = wave_func.real();
            h_psi_old_imag[j * TPB.x * nBlocks.x + i] = wave_func.imag();
            h_psi_new_real[j * TPB.x * nBlocks.x + i] = 0.;
            h_psi_new_imag[j * TPB.x * nBlocks.x + i] = 0.;
            h_potential[j * TPB.x * nBlocks.x + i] = potential_value;
        }
    }

    cudaMemcpy(d_psi_new_real, h_psi_new_real, sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyHostToDevice);
    cudaMemcpy(d_psi_new_imag, h_psi_new_imag, sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyHostToDevice);
    cudaMemcpy(d_psi_old_real, h_psi_old_real, sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyHostToDevice);
    cudaMemcpy(d_psi_old_imag, h_psi_old_imag, sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyHostToDevice);
    cudaMemcpy(d_potential, h_potential, sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyHostToDevice);
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
    cudaMemcpy(h_psi_new_real, d_psi_new_real, sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_psi_new_imag, d_psi_new_imag, sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyDeviceToHost);

    for (int i = 0; i < n_x; ++i)
    {
        for (int j = 0; j < n_y; ++j)
        {
            this->domain->at(i, j, k + 1)->wave_function = h_psi_new_real[j * TPB.x * nBlocks.x + i] + std::complex<double>{0, 1.} * h_psi_new_imag[j * TPB.x * nBlocks.x + i];
        }
    }
}

void FERectPSolver::solve(std::string dir_name="")
{
    int time_length = this->domain->get_num_times();

    for (int k = 0; k < time_length - 1; ++k)
    {
        // std::cout << "Time step: " << k << std::endl;
        this->solve_single_time(k);
        this->domain->normalize(k + 1);
    }
    this->domain->generate_txt_file(std::string{"Forward_Euler_Result"}+dir_name);
}
