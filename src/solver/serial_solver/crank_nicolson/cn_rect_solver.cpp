/**
 * @file cn_rect_solver.cpp
 * @author Minyoung Kim, Gyeonghun Kim
 * @brief Implementation file for serial crank nicolson solver
 * @version 0.1
 * @date 2022-06-05
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "cn_rect_solver.h"

/**
 * @brief Construct a new CNRectSolver::CNRectSolver object
 * 
 * @param g 
 * @param domain 
 */
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
    this->string_info = std::string{"Crank_Nicolson_serial_"};
};

/**
 * @brief initialize guess with forward euler method
 * 
 * @param k 
 */
void CNRectSolver::initialize_guess_with_forward_euler(int k)
{
    delete this->guess;
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

/**
 * @brief Update guess of the spatial(i, j), temporal k point
 * 
 * @param i 
 * @param j 
 * @param k 
 */
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

/**
 * @brief Calculate L2 error between the ccurrent and previous guess
 * 
 * @param k 
 * @return float 
 */
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

/**
 * @brief Solve spatial domain for temporal index k
 * 
 * @param k 
 * @param tolerance 
 * @param max_iter 
 */
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
        std::cout << "At time "<<k <<" Converged failed with error = " << error << std::endl;
    }
}

/**
 * @brief solve equation with crank nicolson method with given max_iter and tolerance.
 * 
 * @param tolerance 
 * @param max_iter 
 * @param dir_name 
 * @param print_info 
 * @param save_data 
 */
void CNRectSolver::solve(float tolerance, int max_iter, std::string dir_name,bool print_info, bool save_data)
{
    int time_length = this->domain->get_num_times();
    //Save initial condition at time = start_time
    if(save_data){
        this->domain->generate_directory_name(this->string_info+dir_name, print_info);
        //Save initial condition
        this->domain->generate_single_txt_file(std::string("Solution_") + std::to_string(0));
    }else{
        this -> domain->update_time();
    }
    
    for (int k = 0; k < time_length-1; ++k)
    {
        //Update kth grid using k-1 th grid
        // std::cout << "time step " << k << std::endl;

        this->initialize_guess_with_forward_euler(k+1);
        this->solve_single_time(k+1, tolerance, max_iter);
        this->domain->normalize(k+1);
        if(save_data){
            this->domain->generate_single_txt_file(std::string("Solution_") + std::to_string(k+1));
        }
        else{
            this->domain->update_time();
        }
    }
    if(print_info){
        this->domain->print_directory_info();
    }
}