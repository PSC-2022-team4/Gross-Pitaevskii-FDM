#include "fe_rect_solver.h"
#include <iostream>
#include <string>
FERectSolver::FERectSolver(
    float g_,
    RectangularDomain *domain_)
    : BaseSolver(g_)
{
    this->domain = domain_;
    this->string_info = std::string{"Forward_Euler_Result_"};
   
};
float FERectSolver::get_potential_value(int i, int j)
{

    return this->domain->potential_grid->at(i, j)->value.real();
}
/**
 * @brief Time differential of phi 
 * 
 * @param i index for x 
 * @param j index for y 
 * @param k index for time(t)
 * @return std::complex<float> time differential at x, y, t 
 */
std::complex<float> FERectSolver::temporal_equation(int i, int j, int k)
{
    //Use five stencil method
    auto point_data = this->domain->at(i, j, k);

    //l,r,d,u denotes left, right, down, up value
    //Check boundary
    auto point_data_l = this->domain->at(i - 1, j, k);
    if (i <= 0)
        point_data_l = this->domain->get_null_gridpt();
    auto point_data_d = this->domain->at(i, j - 1, k);
    if (j <= 0)
        point_data_d = this->domain->get_null_gridpt();
    auto point_data_r = this->domain->at(i + 1, j, k);
    if (i >= (this->domain->get_num_grid_1()) - 1)
        point_data_r = this->domain->get_null_gridpt();
    auto point_data_u = this->domain->at(i, j + 1, k);
    if (j >= (this->domain->get_num_grid_2()) - 1)
        point_data_u = this->domain->get_null_gridpt();

    //potential at x, y
    float V_ij = this->get_potential_value(i, j);
    //this->potential_func(point_data->x, point_data->y);

    //g * |psi(x,y)|^2
    float additional_term = (this->g) * (std::abs(point_data->value)) * (std::abs(point_data->value));

    //Set infinitesimal value
    float dx = this->domain->get_infinitesimal_distance1();
    float dy = this->domain->get_infinitesimal_distance2();
    //df denote time differential of dt (d(psi)/dt)
    // = (laplace - V-g|psi|^2) psi
    std::complex<float> df =
        +((point_data_r->value) + (point_data_l->value) - (point_data->value) * std::complex<float>{2}) / (std::complex<float>{dx * dx}) + ((point_data_u->value) + (point_data_d->value) - (point_data->value) * std::complex<float>{2}) / (std::complex<float>{dy * dy}) - (V_ij + additional_term) * (point_data->value);
    df *= std::complex<float>{0, 1};

    return df;
};

/**
 * @brief Update k+1 th  wave function value with k th values 
 * 
 * @param k 
 */
void FERectSolver::solve_single_time(int k)
{
    int Nx = this->domain->get_num_grid_1();
    int Ny = this->domain->get_num_grid_2();
    float dt = this->domain->get_dt();
    for (int i = 0; i < Nx; ++i)
    {
        for (int j = 0; j < Ny; ++j)
        {
            (this->domain->at(i, j, k + 1)->value) = this->domain->at(i, j, k)->value + dt * this->temporal_equation(i, j, k);
        }
    }
}
void FERectSolver::solve(std::string dir_name, bool print_info, bool save_data)
{
    int time_length = this->domain->get_num_times();
    if(save_data){
        this->domain->generate_directory_name(this->string_info+dir_name, print_info);
        //Save initial condition
        this->domain->generate_single_txt_file(std::string("Solution_") + std::to_string(0));
    }else{
        this -> domain->update_time();
    }
    
    for (int k = 0; k < time_length - 1; ++k)
    {
        this->solve_single_time(k);
        this->domain->normalize(k + 1);
        if(save_data){
            this->domain->generate_single_txt_file(std::string("Solution_") + std::to_string(k + 1));
        }
    }
    if(print_info){
        this->domain->print_directory_info();
    }
}
