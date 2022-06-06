#include "src/serial_solver/crank_nicolson/cn_rect_solver_copy.h"
#include <iostream>
#include <cmath>

CNRectSolver::CNRectSolver(
    std::function<double(double, double)> potential,
    double g,

    RectangularDomain *domain)
    : FERectSolver(potential, g, domain)
{

    this->fe_solver = new FERectSolver(potential, g, domain);
    this->string_info = std::string{"crank_nicolson_serial"};
};

void CNRectSolver::initialize_guess_with_forward_euler(int k)
{
    this->fe_solver->solve_single_time(k-1);
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
            this->guess->at(i, j)->wave_function = this->domain->at(i, j, k)->wave_function;
        }
    }
}

void CNRectSolver::update_guess(int i, int j, int k){
    for (auto i = 0; i < this->domain->get_num_grid_1(); ++i)
    {
        for (auto j = 0; j < this->domain->get_num_grid_2(); ++j)
        {
            this->guess->at(i, j)->wave_function = (0.99 * this->guess->at(i, j)->wave_function +
                                                    0.01 * (this->domain->at(i, j, k - 1)->wave_function + std::complex<double>{0.5, 0.} * this->domain->get_dt() * (this->temporal_equation(i, j, k - 1) + this->temporal_equation(i, j, k)))); 
        }
    }
}

double CNRectSolver::calculate_error(int k){
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

void CNRectSolver::solve_single_time(int k, double tolerance, int max_iter)
{
    double error = 1.;
    bool converged = false;
    int converged_step = 0;
    for (auto iter = 0; iter < max_iter; ++iter)
    {
        //Check convergence 
        if(error < tolerance){
            converged = true;
            converged_step = iter - 1;
            break;
        }
        
        //For each point, update wave function 
        for (auto i = 0; i < this->domain->get_num_grid_1(); ++i)
        {
            for (auto j = 0; j < this->domain->get_num_grid_2(); ++j)
            {   
                this->domain->at(i, j, k)->wave_function = this-> guess->at(i, j)->wave_function;
            }
        }
        
        //for each point, update predicted value
        for (auto i = 0; i < this->domain->get_num_grid_1(); ++i){
            for (auto j = 0; j < this->domain->get_num_grid_2(); ++j){
                update_guess(i, j, k);
                
            }
        }

        error = this->calculate_error(k);
    }
    if(!converged){
        std::cout << "Converged failed with error = " << error << std::endl;
    }
}
void CNRectSolver::solve(double tolerance, int max_iter, std::string dir_name ){

    int time_length = this->domain->get_num_times();
    this->domain ->generate_directory_name(this->string_info);
    for(int k=0; k<time_length-1; ++k){
        //Update kth grid using k-1 th grid 
        std::cout << "time step " << k << std::endl;
        this->initialize_guess_with_forward_euler(k);
        this->solve_single_time(k, tolerance, max_iter);
        this->domain->normalize(k);
        this -> domain->generate_single_txt_file( std::string("probability_") + std::to_string(k));
    }
    
    this->domain->print_directory_info();
}