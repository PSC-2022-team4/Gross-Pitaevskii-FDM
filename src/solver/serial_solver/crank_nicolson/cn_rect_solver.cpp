#include "cn_rect_solver.h"
#include <iostream>
#include <cmath>

CNRectSolver::CNRectSolver(
    // std::function<float(float, float)> potential,
    float g,
    RectangularDomain *domain)
    : BaseSolver(g), domain(domain) //
{
    // this->generate_potential_grid();
    this->fe_solver = new FERectSolver(g, domain); // potential,
    this->string_info = std::string{"crank_nicolson_serial"};
};

// void CNRectSolver::generate_potential_grid()
// {
//     int num_grid_1 = this->domain->get_num_grid_1();
//     int num_grid_2 = this->domain->get_num_grid_2();
//     float x_start = this->domain->at(0,0,0)->x;
//     float y_start = this->domain->at(0,0,0)->y;
//     float x_end = this->domain->at(num_grid_1-1,num_grid_2-1,0)->x;
//     float y_end = this->domain->at(num_grid_1-1,num_grid_2-1,0)->y;
//     this ->potential_grid = RectangularSpatialGrid(num_grid_1, num_grid_2, x_start, x_end, y_start, y_end);
//     for(auto i=0; i<num_grid_1; ++i){
//         for(auto j=0;j<num_grid_2; ++j){
//             auto point = potential_grid.at(i, j);
//             point->value = {this->potential_func(point->x, point->y), 0};
//         }
//     }
// };

/**
 * @brief Time differential of phi 
 * 
 * @param i index for x 
 * @param j index for y 
 * @param k index for time 
 * @return std::complex<float> 
 */
std::complex<float> CNRectSolver::temporal_equation(int i, int j, int k)
{
    auto infinitesimal_distance_1 = this->domain->get_infinitesimal_distance1();
    auto infinitesimal_distance_2 = this->domain->get_infinitesimal_distance2();
    auto point_data = this->domain->at(i, j, k);
    auto point_data_left = this->domain->at(i - 1, j, k);
    auto point_data_right = this->domain->at(i + 1, j, k);
    if (i == 0)
    {
        point_data_left = new GridPoint(0, 0, std::complex<float>{0.});
    }
    if (i == (this->domain->get_num_grid_1() - 1))
    {
        point_data_right = new GridPoint(0, 0, std::complex<float>{0.});
    }
    auto point_data_up = this->domain->at(i, j + 1, k);
    auto point_data_down = this->domain->at(i, j - 1, k);
    if (j == 0)
    {
        point_data_down = new GridPoint(0, 0, std::complex<float>{0.});
    }
    if (j == (this->domain->get_num_grid_2() - 1))
    {
        point_data_up = new GridPoint(0, 0, std::complex<float>{0.});
    }

    auto potential_value = this->domain->potential_grid.at(i, j)->value.real();
    auto laplacian_x = (-2 / (infinitesimal_distance_1 * infinitesimal_distance_1) * point_data->value +
                        1 / (infinitesimal_distance_1 * infinitesimal_distance_1) * point_data_left->value +
                        1 / (infinitesimal_distance_1 * infinitesimal_distance_1) * point_data_right->value);

    auto laplacian_y = (-2 / (infinitesimal_distance_2 * infinitesimal_distance_2) * point_data->value +
                        1 / (infinitesimal_distance_2 * infinitesimal_distance_2) * point_data_down->value +
                        1 / (infinitesimal_distance_2 * infinitesimal_distance_2) * point_data_up->value);

    auto wave_function_abs_square = (point_data->value.real() * point_data->value.real() +
                                     point_data->value.imag() * point_data->value.imag());

    return laplacian_x + laplacian_y - potential_value * point_data->value - this->g * wave_function_abs_square * point_data->value;
}

/**
 * @brief 
 * 
 * @param i 
 * @param j 
 * @return std::complex<float> 
 */
std::complex<float> CNRectSolver::temporal_equation_from_guess(int i, int j)
{
    auto infinitesimal_distance_1 = this->domain->get_infinitesimal_distance1();
    auto infinitesimal_distance_2 = this->domain->get_infinitesimal_distance2();
    auto point_data = this->guess->at(i, j);
    auto point_data_left = this->guess->at(i - 1, j);
    auto point_data_right = this->guess->at(i + 1, j);
    if (i == 0)
    {
        point_data_left = new GridPoint(0, 0, std::complex<float>{0.});
    }
    if (i == (this->domain->get_num_grid_1() - 1))
    {
        point_data_right = new GridPoint(0, 0, std::complex<float>{0.});
    }
    auto point_data_up = this->guess->at(i, j + 1);
    auto point_data_down = this->guess->at(i, j - 1);
    if (j == 0)
    {
        point_data_down = new GridPoint(0, 0, std::complex<float>{0.});
    }
    if (j == (this->domain->get_num_grid_2() - 1))
    {
        point_data_up = new GridPoint(0, 0, std::complex<float>{0.});
    }

    auto potential_value = this->domain->potential_grid.at(i, j)->value.real();
    auto laplacian_x = (-2 / (infinitesimal_distance_1 * infinitesimal_distance_1) * point_data->value +
                        1 / (infinitesimal_distance_1 * infinitesimal_distance_1) * point_data_left->value +
                        1 / (infinitesimal_distance_1 * infinitesimal_distance_1) * point_data_right->value);

    auto laplacian_y = (-2 / (infinitesimal_distance_2 * infinitesimal_distance_2) * point_data->value +
                        1 / (infinitesimal_distance_2 * infinitesimal_distance_2) * point_data_down->value +
                        1 / (infinitesimal_distance_2 * infinitesimal_distance_2) * point_data_up->value);

    auto wave_function_abs_square = (point_data->value.real() * point_data->value.real() +
                                     point_data->value.imag() * point_data->value.imag());
    return laplacian_x + laplacian_y - potential_value * point_data->value - this->g * wave_function_abs_square * point_data->value;
}
void CNRectSolver::initialize_guess_with_forward_euler(int k)
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
            guess->at(i, j)->value = this->domain->at(i, j, k)->value;
        }
    }
}

void CNRectSolver::update_guess(int i, int j, int k)
{
    for (auto i = 0; i < this->domain->get_num_grid_1(); ++i)
    {
        for (auto j = 0; j < this->domain->get_num_grid_2(); ++j)
        {
            this->guess->at(i, j)->value = (0.99f * this->guess->at(i, j)->value +
                                            0.01f * (this->domain->at(i, j, k - 1)->value + std::complex<float>{0, 0.5} * this->domain->get_dt() * (this->temporal_equation(i, j, k - 1) + this->temporal_equation(i, j, k)))); // this->temporal_equation_from_guess(i, j))));
        }
    }
}

float CNRectSolver::calculate_error(int k)
{
    float error = 0.;
    for (auto i = 0; i < this->domain->get_num_grid_1(); ++i)
    {
        for (auto j = 0; j < this->domain->get_num_grid_2(); ++j)
        {
            error += std::pow(std::abs((this->guess->at(i, j)->value - this->domain->at(i, j, k)->value)), 2);
        }
    }
    return error;
}

void CNRectSolver::solve_single_time(int k, float tolerance, int max_iter)
{
    float error = 1.;
    bool converged = false;
    int converged_step = 0;
    for (auto iter = 0; iter < max_iter; ++iter)
    {
        if (error < tolerance)
        {
            converged = true;
            break;
            converged_step = iter - 1;
        }

        for (auto i = 0; i < this->domain->get_num_grid_1(); ++i)
        {
            for (auto j = 0; j < this->domain->get_num_grid_2(); ++j)
            {
                this->domain->at(i, j, k)->value = guess->at(i, j)->value;
            }
        }

        for (auto i = 0; i < this->domain->get_num_grid_1(); ++i)
        {
            for (auto j = 0; j < this->domain->get_num_grid_2(); ++j)
            {
                update_guess(i, j, k);
            }
        }

        error = this->calculate_error(k);
    }
    if (!converged)
    {
        std::cout << "Converged failed with error = " << error << std::endl;
    }
}
void CNRectSolver::solve(float tolerance, int max_iter)
{
    for (auto k = 1; k < this->domain->get_num_times(); ++k)
    {
        // std::cout << "time step " << k << std::endl;
        this->initialize_guess_with_forward_euler(k);
        this->solve_single_time(k, tolerance, max_iter);
    }
    // this->domain->generate_txt_file(string_info);
}