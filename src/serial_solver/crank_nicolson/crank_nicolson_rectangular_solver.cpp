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
    auto point_data_right = this->domain->at(i + 1, j, k);
    if (i == 0){
        point_data_left = new GridPoint(0, 0, std::complex<double>{0.});
    }
    if (i == (this->domain->get_num_grid_1() - 1))
    {
        point_data_right = new GridPoint(0, 0, std::complex<double>{0.});
    }
    auto point_data_up = this->domain->at(i, j+1, k);
    auto point_data_down = this->domain->at(i, j-1, k);
    if (j == 0)
    {
        point_data_down = new GridPoint(0, 0, std::complex<double>{0.});
    }
    if (j == (this->domain->get_num_grid_2() - 1))
    {
        point_data_up = new GridPoint(0, 0, std::complex<double>{0.});
    }

    auto potential_value = this->potential_func(point_data->x, point_data->y);
    auto laplacian_x = (-2. / (infinitesimal_distance_1 * infinitesimal_distance_1) * point_data->wave_function +
                        1. / (infinitesimal_distance_1 * infinitesimal_distance_1) * point_data_left->wave_function +
                        1. / (infinitesimal_distance_1 * infinitesimal_distance_1) * point_data_right->wave_function);


    auto laplacian_y = (-2. / (infinitesimal_distance_2 * infinitesimal_distance_2) * point_data->wave_function + 
                        1. / (infinitesimal_distance_2 * infinitesimal_distance_2) * point_data_down->wave_function + 
                        1. / (infinitesimal_distance_2 * infinitesimal_distance_2) * point_data_up->wave_function);

    auto wave_function_abs_square = (
        point_data->wave_function.real() * point_data->wave_function.real() + 
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
std::complex<double> CrankNicolsonRectangularSolver::temporal_equation_from_guess(int i, int j)
{
    auto infinitesimal_distance_1 = this->domain->get_infinitesimal_distance1();
    auto infinitesimal_distance_2 = this->domain->get_infinitesimal_distance2();
    auto point_data = this->old_guess->at(i, j);
    auto point_data_left = this->old_guess->at(i - 1, j);
    auto point_data_right = this->old_guess->at(i + 1, j);
    if (i == 0)
    {
        point_data_left = new GridPoint(0, 0, std::complex<double>{0.});
    }
    if (i == (this->domain->get_num_grid_1() - 1))
    {
        point_data_right = new GridPoint(0, 0, std::complex<double>{0.});
    }
    auto point_data_up = this->old_guess->at(i, j - 1);
    auto point_data_down = this->old_guess->at(i, j + 1);
    if (j == 0)
    {
        point_data_down = new GridPoint(0, 0, std::complex<double>{0.});
    }
    if (j == (this->domain->get_num_grid_2() - 1))
    {
        point_data_up = new GridPoint(0, 0, std::complex<double>{0.});
    }

    auto potential_value = this->potential_func(point_data->x, point_data->y);
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
void CrankNicolsonRectangularSolver::initialize_guess_with_forward_euler(int k)
{
    // Update old_guess with Forward Euler result
}

void CrankNicolsonRectangularSolver::update_guess(int i, int j, int k){

}
void CrankNicolsonRectangularSolver::solve_single_time(int k, double tolerance, int max_iter){
    double error = 1.;
    for (auto iter = 0; iter < max_iter; ++iter)
    {
        if(error < tolerance)
            break;
                

    }
}
void CrankNicolsonRectangularSolver::solve(double tolerance, int max_iter){
    this->applyInitialCondition();
    for (auto k = 1; k < this->domain->get_num_times(); ++k)
    {
        this->initialize_guess_with_forward_euler(k);
        this->solve_single_time(k, tolerance, max_iter);
    }
}