#include "cn_rect_solver.h"
#include <iostream>
#include <cmath>

CNRectSolver::CNRectSolver(
    float g,
    RectangularDomain *domain)
    : FERectSolver(g, domain) //
{
    // this->generate_potential_grid();
    this->guess = new RectangularSpatialGrid(
        this->domain->get_num_grid_1(),
        this->domain->get_num_grid_2(),
        this->domain->get_x_start(),
        this->domain->get_x_end(),
        this->domain->get_y_start(),
        this->domain->get_y_end());
    this->string_info = std::string{"crank_nicolson_serial"};
};

void CNRectSolver::initialize_guess_with_forward_euler(int k)
{
    free(this->guess);
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

            //At first, initial guess is the results of fe algorithm
            this->domain->at(i, j, k)->value = this->domain->at(i, j, k - 1)->value + this->domain->get_dt() * (this->temporal_equation(i, j, k - 1));

            this->domain->normalize(k);
            this->guess->at(i, j)->value = this->domain->at(i, j, k)->value;
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
                                            0.01f * (this->domain->at(i, j, k - 1)->value + std::complex<float>{0.5, 0.} * this->domain->get_dt() * (this->temporal_equation(i, j, k - 1) + this->temporal_equation(i, j, k)))); // this->temporal_equation_from_guess(i, j))));
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
            converged_step = iter - 1;
            break;
        }
        //For each point, update wave function
        for (auto i = 0; i < this->domain->get_num_grid_1(); ++i)
        {
            for (auto j = 0; j < this->domain->get_num_grid_2(); ++j)
            {
                this->domain->at(i, j, k)->value = this->guess->at(i, j)->value;
            }
        }

        //for each point, update predicted value
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
void CNRectSolver::solve(float tolerance, int max_iter, std::string dir_name)
{
    int time_length = this->domain->get_num_times();
    //Save initial condition at time = start_time
    this->domain->generate_directory_name(this->string_info + dir_name);
    this->domain->generate_single_txt_file(std::string("probability_") + std::to_string(0));

    for (int k = 1; k < time_length; ++k)
    {
        //Update kth grid using k-1 th grid
        // std::cout << "time step " << k << std::endl;
        this->initialize_guess_with_forward_euler(k);
        this->solve_single_time(k, tolerance, max_iter);
        this->domain->normalize(k);
        this->domain->generate_single_txt_file(std::string("probability_") + std::to_string(k));
    }

    this->domain->print_directory_info();
}