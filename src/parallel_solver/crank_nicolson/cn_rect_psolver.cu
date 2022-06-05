#include "src/parallel_solver/crank_nicolson/cn_rect_psolver.cuh"
#include <iostream>
#include <cmath>

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
                                 double tau)
{
    __shared__ double tile_old_real[nTx][nTy];
    __shared__ double tile_old_imag[nTx][nTy];
    __shared__ double tile_new_real_trial[nTx][nTy];
    __shared__ double tile_new_imag_trial[nTx][nTy];
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
    tile_new_real_trial[thread_x][thread_y] = psi_new_real_trial[j * striding + i];
    tile_new_imag_trial[thread_x][thread_y] = psi_new_imag_trial[j * striding + i];
    tile_potential[thread_x][thread_y] = potential[j * striding + i];

    __syncthreads();

    //  Update tile position
    psi_new_real[j * striding + i] += 0.5 * ((1 - tile_potential[thread_x][thread_y] - g * (tile_old_real[thread_x][thread_y] * tile_old_real[thread_x][thread_y] + tile_old_imag[thread_x][thread_y] * tile_old_imag[thread_x][thread_y])) * tile_old_real[thread_x][thread_y] + (2 * tau / (h_x * h_x) + 2 * tau / (h_y * h_y)) * tile_old_imag[thread_x][thread_y]) *
                                      (i >= 0) * (i < n_x) * (j >= 0) * (j < n_y);
    psi_new_imag[j * striding + i] += 0.5 * ((1 - tile_potential[thread_x][thread_y] - g * (tile_old_real[thread_x][thread_y] * tile_old_real[thread_x][thread_y] + tile_old_imag[thread_x][thread_y] * tile_old_imag[thread_x][thread_y])) * tile_old_imag[thread_x][thread_y] - (2 * tau / (h_x * h_x) + 2 * tau / (h_y * h_y)) * tile_old_real[thread_x][thread_y]) *
                                      (i >= 0) * (i < n_x) * (j >= 0) * (j < n_y);
    psi_new_real[j * striding + i] += 0.5 * ((1 - tile_potential[thread_x][thread_y] - g * (tile_new_real_trial[thread_x][thread_y] * tile_new_real_trial[thread_x][thread_y] + tile_new_imag_trial[thread_x][thread_y] * tile_new_imag_trial[thread_x][thread_y])) * tile_new_real_trial[thread_x][thread_y] + (2 * tau / (h_x * h_x) + 2 * tau / (h_y * h_y)) * tile_new_imag_trial[thread_x][thread_y]) *
                                      (i >= 0) * (i < n_x) * (j >= 0) * (j < n_y);
    psi_new_imag[j * striding + i] += 0.5 * ((1 - tile_potential[thread_x][thread_y] - g * (tile_new_real_trial[thread_x][thread_y] * tile_new_real_trial[thread_x][thread_y] + tile_new_imag_trial[thread_x][thread_y] * tile_new_imag_trial[thread_x][thread_y])) * tile_new_imag_trial[thread_x][thread_y] - (2 * tau / (h_x * h_x) + 2 * tau / (h_y * h_y)) * tile_new_real_trial[thread_x][thread_y]) *
                                      (i >= 0) * (i < n_x) * (j >= 0) * (j < n_y);
    // Update left position
    psi_new_real[j * striding + i - 1] += 0.5 * (-(tau / (h_x * h_x)) * tile_old_imag[thread_x][thread_y]) * (i > 0) * (i < n_x) * (j >= 0) * (j < n_y);
    psi_new_imag[j * striding + i - 1] += 0.5 * ((tau / (h_x * h_x)) * tile_old_real[thread_x][thread_y]) * (i > 0) * (i < n_x) * (j >= 0) * (j < n_y);
    psi_new_real[j * striding + i - 1] += 0.5 * (-(tau / (h_x * h_x)) * tile_new_imag_trial[thread_x][thread_y]) * (i > 0) * (i < n_x) * (j >= 0) * (j < n_y);
    psi_new_imag[j * striding + i - 1] += 0.5 * ((tau / (h_x * h_x)) * tile_new_real_trial[thread_x][thread_y]) * (i > 0) * (i < n_x) * (j >= 0) * (j < n_y);

    // Update right position
    psi_new_real[j * striding + i + 1] += 0.5 * (-(tau / (h_x * h_x)) * tile_old_imag[thread_x][thread_y]) * (i >= 0) * (i < (n_x - 1)) * (j >= 0) * (j < n_y);
    psi_new_imag[j * striding + i + 1] += 0.5 * ((tau / (h_x * h_x)) * tile_old_real[thread_x][thread_y]) * (i >= 0) * (i < (n_x - 1)) * (j >= 0) * (j < n_y);
    psi_new_real[j * striding + i + 1] += 0.5 * (-(tau / (h_x * h_x)) * tile_new_imag_trial[thread_x][thread_y]) * (i >= 0) * (i < (n_x - 1)) * (j >= 0) * (j < n_y);
    psi_new_imag[j * striding + i + 1] += 0.5 * ((tau / (h_x * h_x)) * tile_new_real_trial[thread_x][thread_y]) * (i >= 0) * (i < (n_x - 1)) * (j >= 0) * (j < n_y);

    // Update down position
    psi_new_real[(j - 1) * striding + i] += 0.5 * (-(tau / (h_y * h_y)) * tile_old_imag[thread_x][thread_y]) * (i >= 0) * (i < n_x) * (j > 0) * (j < n_y);
    psi_new_imag[(j - 1) * striding + i] += 0.5 * ((tau / (h_y * h_y)) * tile_old_real[thread_x][thread_y]) * (i >= 0) * (i < n_x) * (j > 0) * (j < n_y);
    psi_new_real[(j - 1) * striding + i] += 0.5 * (-(tau / (h_y * h_y)) * tile_new_imag_trial[thread_x][thread_y]) * (i >= 0) * (i < n_x) * (j > 0) * (j < n_y);
    psi_new_imag[(j - 1) * striding + i] += 0.5 * ((tau / (h_y * h_y)) * tile_new_real_trial[thread_x][thread_y]) * (i >= 0) * (i < n_x) * (j > 0) * (j < n_y);

    // Update up position
    psi_new_real[(j + 1) * striding + i] += 0.5 * (-(tau / (h_y * h_y)) * tile_old_imag[thread_x][thread_y]) * (i >= 0) * (i < n_x) * (j >= 0) * (j < (n_y - 1));
    psi_new_imag[(j + 1) * striding + i] += 0.5 * ((tau / (h_y * h_y)) * tile_old_real[thread_x][thread_y]) * (i >= 0) * (i < n_x) * (j >= 0) * (j < (n_y - 1));
    psi_new_real[(j + 1) * striding + i] += 0.5 * (-(tau / (h_y * h_y)) * tile_new_imag_trial[thread_x][thread_y]) * (i >= 0) * (i < n_x) * (j >= 0) * (j < (n_y - 1));
    psi_new_imag[(j + 1) * striding + i] += 0.5 * ((tau / (h_y * h_y)) * tile_new_real_trial[thread_x][thread_y]) * (i >= 0) * (i < n_x) * (j >= 0) * (j < (n_y - 1));
}

__global__ void calculate_error_on_device(double *psi_1_real,
                                          double *psi_1_imag,
                                          double *psi_2_real,
                                          double *psi_2_imag,
                                          double *error)
{
    __shared__ double tile_psi_1_real[nTx][nTy];
    __shared__ double tile_psi_1_imag[nTx][nTy];
    __shared__ double tile_psi_2_real[nTx][nTy];
    __shared__ double tile_psi_2_imag[nTx][nTy];

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
    *error += (tile_psi_1_real[thread_x][thread_y] - tile_psi_2_real[thread_x][thread_y]) * (tile_psi_1_real[thread_x][thread_y] - tile_psi_2_real[thread_x][thread_y]);
    *error += (tile_psi_1_imag[thread_x][thread_y] - tile_psi_2_imag[thread_x][thread_y]) * (tile_psi_1_imag[thread_x][thread_y] - tile_psi_2_imag[thread_x][thread_y]);
    printf("%f", *error);
}

CNRectPSolver::CNRectPSolver(
    std::function<double(double, double)> potential,
    double g,

    RectangularDomain *domain)
    : BaseSolver(potential, g), domain(domain)
{
    this->generate_potential_grid();
    this->fe_solver = new FERectSolver(potential, g, domain);
};

void CNRectPSolver::generate_potential_grid()
{
    int num_grid_1 = this->domain->get_num_grid_1();
    int num_grid_2 = this->domain->get_num_grid_2();
    double x_start = this->domain->at(0, 0, 0)->x;
    double y_start = this->domain->at(0, 0, 0)->y;
    double x_end = this->domain->at(num_grid_1 - 1, num_grid_2 - 1, 0)->x;
    double y_end = this->domain->at(num_grid_1 - 1, num_grid_2 - 1, 0)->y;
    this->potential_grid = RectangularSpatialGrid(num_grid_1, num_grid_2, x_start, x_end, y_start, y_end);
    for (auto i = 0; i < num_grid_1 - 1; ++i)
    {
        for (auto j = 0; j < num_grid_2; ++j)
        {
            auto point = potential_grid.at(i, j);
            point->wave_function = {this->potential_func(point->x, point->y), 0};
        }
    }
};

/**
 * @brief Time differential of phi
 *
 * @param i index for x
 * @param j index for y
 * @param k index for time
 * @return std::complex<double>
 */
std::complex<double> CNRectPSolver::temporal_equation(int i, int j, int k)
{
    auto infinitesimal_distance_1 = this->domain->get_infinitesimal_distance1();
    auto infinitesimal_distance_2 = this->domain->get_infinitesimal_distance2();
    auto point_data = this->domain->at(i, j, k);
    auto point_data_left = this->domain->at(i - 1, j, k);
    auto point_data_right = this->domain->at(i + 1, j, k);
    if (i == 0)
    {
        point_data_left = new GridPoint(0, 0, std::complex<double>{0.});
    }
    if (i == (this->domain->get_num_grid_1() - 1))
    {
        point_data_right = new GridPoint(0, 0, std::complex<double>{0.});
    }
    auto point_data_up = this->domain->at(i, j + 1, k);
    auto point_data_down = this->domain->at(i, j - 1, k);
    if (j == 0)
    {
        point_data_down = new GridPoint(0, 0, std::complex<double>{0.});
    }
    if (j == (this->domain->get_num_grid_2() - 1))
    {
        point_data_up = new GridPoint(0, 0, std::complex<double>{0.});
    }

    auto potential_value = this->potential_grid.at(i, j)->wave_function.real();
    auto laplacian_x = (-2. / (infinitesimal_distance_1 * infinitesimal_distance_1) * point_data->wave_function +
                        1. / (infinitesimal_distance_1 * infinitesimal_distance_1) * point_data_left->wave_function +
                        1. / (infinitesimal_distance_1 * infinitesimal_distance_1) * point_data_right->wave_function);

    auto laplacian_y = (-2. / (infinitesimal_distance_2 * infinitesimal_distance_2) * point_data->wave_function +
                        1. / (infinitesimal_distance_2 * infinitesimal_distance_2) * point_data_down->wave_function +
                        1. / (infinitesimal_distance_2 * infinitesimal_distance_2) * point_data_up->wave_function);

    auto wave_function_abs_square = (point_data->wave_function.real() * point_data->wave_function.real() +
                                     point_data->wave_function.imag() * point_data->wave_function.imag());

    return laplacian_x + laplacian_y - potential_value * point_data->wave_function - this->g * wave_function_abs_square * point_data->wave_function;
}

/**
 * @brief
 *
 * @param i
 * @param j
 * @return std::complex<double>
 */
std::complex<double> CNRectPSolver::temporal_equation_from_guess(int i, int j)
{
    auto infinitesimal_distance_1 = this->domain->get_infinitesimal_distance1();
    auto infinitesimal_distance_2 = this->domain->get_infinitesimal_distance2();
    auto point_data = this->guess->at(i, j);
    auto point_data_left = this->guess->at(i - 1, j);
    auto point_data_right = this->guess->at(i + 1, j);
    if (i == 0)
    {
        point_data_left = new GridPoint(0, 0, std::complex<double>{0.});
    }
    if (i == (this->domain->get_num_grid_1() - 1))
    {
        point_data_right = new GridPoint(0, 0, std::complex<double>{0.});
    }
    auto point_data_up = this->guess->at(i, j + 1);
    auto point_data_down = this->guess->at(i, j - 1);
    if (j == 0)
    {
        point_data_down = new GridPoint(0, 0, std::complex<double>{0.});
    }
    if (j == (this->domain->get_num_grid_2() - 1))
    {
        point_data_up = new GridPoint(0, 0, std::complex<double>{0.});
    }

    auto potential_value = this->potential_grid.at(i, j)->wave_function.real();
    auto laplacian_x = (-2. / (infinitesimal_distance_1 * infinitesimal_distance_1) * point_data->wave_function +
                        1. / (infinitesimal_distance_1 * infinitesimal_distance_1) * point_data_left->wave_function +
                        1. / (infinitesimal_distance_1 * infinitesimal_distance_1) * point_data_right->wave_function);

    auto laplacian_y = (-2. / (infinitesimal_distance_2 * infinitesimal_distance_2) * point_data->wave_function +
                        1. / (infinitesimal_distance_2 * infinitesimal_distance_2) * point_data_down->wave_function +
                        1. / (infinitesimal_distance_2 * infinitesimal_distance_2) * point_data_up->wave_function);

    auto wave_function_abs_square = (point_data->wave_function.real() * point_data->wave_function.real() +
                                     point_data->wave_function.imag() * point_data->wave_function.imag());
    return laplacian_x + laplacian_y - potential_value * point_data->wave_function - this->g * wave_function_abs_square * point_data->wave_function;
}
void CNRectPSolver::initialize_guess_with_forward_euler(int k)
{
    this->fe_solver->solve_single_time(k - 1);
    this->guess = new RectangularSpatialGrid(
        this->domain->get_num_grid_1(),
        this->domain->get_num_grid_2(),
        this->domain->get_x_start(),
        this->domain->get_x_end(),
        this->domain->get_y_start(),
        this->domain->get_y_end());

    for (auto i = 0; i < this->domain->get_num_grid_1(); ++i)
    {
        for (auto j = 0; j < this->domain->get_num_grid_2(); ++j)
        {
            guess->at(i, j)->wave_function = this->domain->at(i, j, k)->wave_function;
        }
    }
}

void CNRectPSolver::update_guess(int i, int j, int k)
{
    for (auto i = 0; i < this->domain->get_num_grid_1(); ++i)
    {
        for (auto j = 0; j < this->domain->get_num_grid_2(); ++j)
        {
            this->guess->at(i, j)->wave_function = (0.99 * this->guess->at(i, j)->wave_function +
                                                    0.01 * (this->domain->at(i, j, k - 1)->wave_function + std::complex<double>{0, 0.5} * this->domain->get_dt() * (this->temporal_equation(i, j, k - 1) + this->temporal_equation(i, j, k)))); // this->temporal_equation_from_guess(i, j))));
        }
    }
}

double CNRectPSolver::calculate_error(int k)
{
    double error = 0.;
    for (auto i = 0; i < this->domain->get_num_grid_1(); ++i)
    {
        for (auto j = 0; j < this->domain->get_num_grid_2(); ++j)
        {
            error += std::pow(std::abs((this->guess->at(i, j)->wave_function - this->domain->at(i, j, k)->wave_function)), 2);
        }
    }
    return error;
}

void CNRectPSolver::solve_single_time(int k, double tolerance, int max_iter)
{
    int n_x = this->domain->get_num_grid_1();
    int n_y = this->domain->get_num_grid_2();
    double dt = this->domain->get_dt();

    double error = 1.;
    double *d_error;
    bool converged = false;
    int converged_step = 0;

    dim3 TPB(nTx, nTy);
    dim3 nBlocks(n_x / nTx + (n_x % nTx != 0), n_y / nTy + (n_y % nTy != 0));
    double *h_psi_old_real = (double *)malloc(sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    double *h_psi_old_imag = (double *)malloc(sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    double *h_psi_new_real_trial = (double *)malloc(sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    double *h_psi_new_imag_trial = (double *)malloc(sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    double *h_psi_new_real = (double *)malloc(sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    double *h_psi_new_imag = (double *)malloc(sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    double *h_potential = (double *)malloc(sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    double *d_psi_old_real, *d_psi_old_imag, *d_psi_new_real, *d_psi_new_real_trial, *d_psi_new_imag_trial, *d_psi_new_imag, *d_potential;
    cudaMalloc((double **)&d_psi_old_real, sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    cudaMalloc((double **)&d_psi_old_imag, sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    cudaMalloc((double **)&d_psi_new_real_trial, sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    cudaMalloc((double **)&d_psi_new_imag_trial, sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    cudaMalloc((double **)&d_psi_new_real, sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    cudaMalloc((double **)&d_psi_new_imag, sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    cudaMalloc((double **)&d_potential, sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y);
    cudaMalloc((double **)&d_error, sizeof(double));
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
            h_psi_new_real_trial[j * TPB.x * nBlocks.x + i] = wave_func.real();
            h_psi_new_imag_trial[j * TPB.x * nBlocks.x + i] = wave_func.imag();
            h_psi_new_real[j * TPB.x * nBlocks.x + i] = wave_func.real();
            h_psi_new_imag[j * TPB.x * nBlocks.x + i] = wave_func.imag();
            h_potential[j * TPB.x * nBlocks.x + i] = potential_value;
        }
    }

    cudaMemcpy(d_psi_new_real, h_psi_new_real, sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyHostToDevice);
    cudaMemcpy(d_psi_new_imag, h_psi_new_imag, sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyHostToDevice);
    cudaMemcpy(d_psi_new_real_trial, h_psi_new_real_trial, sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyHostToDevice);
    cudaMemcpy(d_psi_new_imag_trial, h_psi_new_imag_trial, sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyHostToDevice);
    cudaMemcpy(d_psi_old_real, h_psi_old_real, sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyHostToDevice);
    cudaMemcpy(d_psi_old_imag, h_psi_old_imag, sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyHostToDevice);
    cudaMemcpy(d_potential, h_potential, sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyHostToDevice);
    cudaMemcpy(d_error, &error, sizeof(double), cudaMemcpyHostToDevice);

    cudaDeviceSynchronize();

    for (auto iter = 0; iter < max_iter; ++iter)
    {
        if (error < tolerance)
        {
            converged = true;
            converged_step = iter - 1;
            break;
        }
        cudaMemcpy(d_psi_new_real_trial, d_psi_new_real, sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyDeviceToDevice);
        cudaMemcpy(d_psi_new_imag_trial, d_psi_new_imag, sizeof(double) * TPB.x * nBlocks.x * TPB.y * nBlocks.y, cudaMemcpyDeviceToDevice);
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
            this->domain->get_dt());
        error = 0.;
        cudaMemcpy(d_error, &error, sizeof(double), cudaMemcpyHostToDevice);
        calculate_error_on_device<<<nBlocks, TPB>>>(d_psi_new_real_trial,
                                                    d_psi_new_imag_trial,
                                                    d_psi_new_real,
                                                    d_psi_new_imag,
                                                    &error);
        cudaMemcpy(&error, d_error, sizeof(double), cudaMemcpyDeviceToHost);
        std::cout << error << std::endl;
    }
    if (!converged)
    {
        std::cout << "Converged failed with error = " << error << std::endl;
    }
}
void CNRectPSolver::solve(double tolerance, int max_iter)
{
    for (auto k = 1; k < this->domain->get_num_times(); ++k)
    {
        std::cout << "time step " << k << std::endl;
        this->initialize_guess_with_forward_euler(k);
        this->solve_single_time(k, tolerance, max_iter);
    }
    this->domain->generate_txt_file(std::string{"Forward_Euler_Result"});
}