#include "src/serial_solver/crank_nicolson/crank_nicolson_rectangular_solver.h"

CrankNicolsonRectangularSolver::CrankNicolsonRectangularSolver(
    InitialCondition initialCondition,
    std::function<double(double, double)> potential,
    double g,
    RectangularDomain rectangularDomain)
    : BaseSolver(initialCondition, potential, g, rectangularDomain)
{
    // Cast basedomain to rectangularDomain
    (RectangularDomain) * (this->domain);
    };

void CrankNicolsonRectangularSolver::applyInitialCondition()
{
    this->initialCondition.assign_to_domain(this->domain);
}

/**
 * @brief Time differential of phi 
 * 
 * @param i index for x 
 * @param j index for y 
 * @param k index for time 
 * @return std::complex<double> 
 */
std::complex<double> CrankNicolsonRectangularSolver::temporal_equation(int i, int j, int k)
{
    auto infinitesimal_distance_1 = this->domain->get_infinitesimal_distance1();
    auto infinitesimal_distance_2 = this->domain->get_infinitesimal_distance2();
    auto point_data = this->domain->at(i, j, k);
    auto point_data_left = this->domain->at(i - 1, j, k);
    auto point_data_right = this->domain->at(i+1, j, k);
    auto point_data_up = this->domain->at(i, j-1, k);
    auto point_data_down = this->domain->at(i, j+1, k);

    auto potential_value = this->potential_func(point_data->x, point_data->y);
    auto laplacian_x = -2. / (infinitesimal_distance_1 * infinitesimal_distance_1) * point_data->wave_function;
    if (i > 0){
        laplacian_x += 1. / (infinitesimal_distance_1 * infinitesimal_distance_1) * point_data_left->wave_function;
    }
    if (i < (this->domain->get_num_grid_1()-1)){
        laplacian_x += 1. / (infinitesimal_distance_1 * infinitesimal_distance_1) * point_data_right->wave_function;
    }

    auto laplacian_y = -2. / (infinitesimal_distance_2 * infinitesimal_distance_2) * point_data->wave_function;
    if (j > 0)
    {
        laplacian_y += 1. / (infinitesimal_distance_2 * infinitesimal_distance_2) * point_data_down->wave_function;
    }
    if (j < (this->domain->get_num_grid_2()-1))
    {
        laplacian_y += 1. / (infinitesimal_distance_2 * infinitesimal_distance_2) * point_data_up->wave_function;
    }
    auto wave_function_abs_square = (point_data->wave_function.real() * point_data->wave_function.real() + point_data->wave_function.imag() * point_data->wave_function.imag());
    return laplacian_x + laplacian_y - potential_value * point_data->wave_function - this->g * wave_function_abs_square * point_data->wave_function;
}

/**
 * @brief 
 * 
 * @param i 
 * @param j 
 * @return std::complex<double> 
 */
std::complex<double> CrankNicolsonRectangularSolver::temporal_equation_from_guess(int i, int j)
{
    auto infinitesimal_distance_1 = this->domain->get_infinitesimal_distance1();
    auto infinitesimal_distance_2 = this->domain->get_infinitesimal_distance2();
    auto point_data = this->old_guess->at(i, j);
    auto point_data_left = this->old_guess->at(i - 1, j);
    auto point_data_right = this->old_guess->at(i + 1, j);
    auto point_data_up = this->old_guess->at(i, j - 1);
    auto point_data_down = this->old_guess->at(i, j + 1);

    auto potential_value = this->potential_func(point_data->x, point_data->y);
    auto laplacian_x = -2. / (infinitesimal_distance_1 * infinitesimal_distance_1) * point_data->wave_function;
    if (i > 0)
    {
        laplacian_x += 1. / (infinitesimal_distance_1 * infinitesimal_distance_1) * point_data_left->wave_function;
    }
    if (i < (this->domain->get_num_grid_1() - 1))
    {
        laplacian_x += 1. / (infinitesimal_distance_1 * infinitesimal_distance_1) * point_data_right->wave_function;
    }

    auto laplacian_y = -2. / (infinitesimal_distance_2 * infinitesimal_distance_2) * point_data->wave_function;
    if (j > 0)
    {
        laplacian_y += 1. / (infinitesimal_distance_2 * infinitesimal_distance_2) * point_data_down->wave_function;
    }
    if (j < (this->domain->get_num_grid_2() - 1))
    {
        laplacian_y += 1. / (infinitesimal_distance_2 * infinitesimal_distance_2) * point_data_up->wave_function;
    }
    auto wave_function_abs_square = (point_data->wave_function.real() * point_data->wave_function.real() + point_data->wave_function.imag() * point_data->wave_function.imag());
    return laplacian_x + laplacian_y - potential_value * point_data->wave_function - this->g * wave_function_abs_square * point_data->wave_function;
}

