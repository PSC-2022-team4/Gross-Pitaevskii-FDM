#include "fe_rect_solver.h"
#include <iostream>
#include <string>
FERectSolver::FERectSolver(
    std::function<float(float, float)> potential_, 
    float g_, 
    RectangularDomain* domain_)
    :BaseSolver(potential_, g_){
        this-> domain = domain_;
        this->generate_potential_grid();
        this->string_info = std::string{"Forward_Euler_Result"};
        string_info += "_"+std::to_string(g_);
};
void FERectSolver::generate_potential_grid(){
    int num_grid_1 = this->domain->get_num_grid_1();
    int num_grid_2 = this->domain->get_num_grid_2();
    float x_start = this->domain->at(0,0,0)->x;
    float y_start = this->domain->at(0,0,0)->y;
    float x_end = this->domain->at(num_grid_1-1,num_grid_2-1,0)->x;
    float y_end = this->domain->at(num_grid_1-1,num_grid_2-1,0)->y;
    this ->potential_grid = RectangularSpatialGrid(num_grid_1, num_grid_2, x_start, x_end, y_start, y_end);
    for(auto i=0; i<num_grid_1; ++i){
        for(auto j=0;j<num_grid_2; ++j){
            auto point = potential_grid.at(i, j);
            point->wave_function = {this->potential_func(point->x, point->y), 0};
        }
    }
};
float FERectSolver::get_potential_value(int i, int j ){
    return this->potential_grid.at(i, j)->wave_function.real();
}
/**
 * @brief Time differential of phi 
 * 
 * @param i index for x 
 * @param j index for y 
 * @param k index for time(t)
 * @return std::complex<float> time differential at x, y, t 
 */
std::complex<float> FERectSolver::temporal_equation(int i, int j, int k){
    //Use five stencil method 
    auto point_data = this-> domain->at(i, j, k);

    //l,r,d,u denotes left, right, down, up value 
    //Check boundary 
    auto point_data_l = this-> domain->at(i-1, j, k);
    if(i <= 0)
        point_data_l =this->domain->get_null_gridpt();
    auto point_data_d = this-> domain->at(i, j-1, k);
    if(j<=0 )
        point_data_d = this->domain->get_null_gridpt();
    auto point_data_r = this-> domain->at(i+1, j, k);
    if(i>= (this->domain->get_num_grid_1())-1)
        point_data_r =this->domain->get_null_gridpt();
    auto point_data_u = this-> domain->at(i, j+1, k);
    if(j >= (this->domain->get_num_grid_2())-1)
        point_data_u = this->domain->get_null_gridpt();
    
    //potential at x, y 
    float V_ij = this->get_potential_value(i,j);
    //this->potential_func(point_data->x, point_data->y);
    
    //g * |psi(x,y)|^2 
    float additional_term = (this-> g) * (std::abs(point_data->wave_function))*(std::abs(point_data->wave_function));
    
    
    //Set infinitesimal value 
    float dx = this->domain->get_infinitesimal_distance1();
    float dy = this->domain->get_infinitesimal_distance2();
    //df denote time differential of dt (d(psi)/dt)
    // = (laplace - V-g|psi|^2) psi
    std::complex<float> df = 
        +((point_data_r->wave_function)+(point_data_l->wave_function)-(point_data->wave_function) * std::complex<float>{2})/ (std::complex<float>{dx*dx})
        +((point_data_u->wave_function)+(point_data_d->wave_function)-(point_data->wave_function) * std::complex<float>{2})/ (std::complex<float>{dy*dy})
        - (V_ij+additional_term) * (point_data->wave_function);
    df *= std::complex<float>{0,1}; 
    

    return df; 
};

/**
 * @brief Update k+1 th  wave function value with k th values 
 * 
 * @param k 
 */
void FERectSolver::solve_single_time(int k)
{   int Nx = this->domain->get_num_grid_1();
    int Ny = this->domain->get_num_grid_2();
    float dt = this->domain->get_dt();
    std::cout<<"At time "<<k << " called to update "<<k+1<<"th grid"<<std::endl;
    //std::cout<<this->domain->at(10, 10, k)->wave_function <<"/dt= "<<this->temporal_equation(10,10,k)<<std::endl;
    for(int i=0; i<Nx; ++i){
        for(int j=0; j<Ny; ++j){
            (this->domain->at(i, j, k+1)->wave_function) = this->domain->at(i, j, k)->wave_function + dt * this->temporal_equation(i,j,k);
        }
    }
}
void FERectSolver::solve(){
    
    int time_length = this->domain->get_num_times();
    this->domain ->generate_directory_name(this->string_info);
    //Save initial condition
    this -> domain->generate_single_txt_file( std::string("probability_") + std::to_string(0));
    
    for(int k=0; k<time_length-1; ++k){
        //Update k+1 th grid using k th grid 
        this->solve_single_time(k);
        this->domain->normalize(k+1);
        this ->domain->generate_single_txt_file( std::string("probability_") + std::to_string(k+1));
    }

    this->domain->print_directory_info();

}


