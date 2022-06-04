#include "forward_euler_rectangular_solver.h"

ForwardEulerRectangularSolver::ForwardEulerRectangularSolver(
    InitialCondition initialCondition_, 
    std::function<double(double, double)> potential_, 
    double g_, 
    RectangularDomain* domain_)
    :BaseSolver(initialCondition_, potential_, g_), domain(domain_)
    {
        
        //Cast basedomain to rectangularDomain
        ;
    };

/**
 * @brief Time differential of phi 
 * 
 * @param i index for x 
 * @param j index for y 
 * @param k index for time(t)
 * @return std::complex<double> time differential at x, y, t 
 */
std::complex<double> ForwardEulerRectangularSolver::temporal_equation(int i, int j, int k){
    //Use five stencil method 
    auto point_data = this-> domain->at(i, j, k);

    //l,r,d,u denotes left, right, down, up value 
    //Check boundary 
    auto point_data_l = this-> domain->at(i-1, j, k);
    if(i <= 0)
        point_data_l = new GridPoint(0., 0., std::complex<double>{0,0});
    auto point_data_d = this-> domain->at(i, j-1, k);
    if(j<=0 )
        point_data_d = new GridPoint(0., 0., std::complex<double>{0,0});
    auto point_data_r = this-> domain->at(i+1, j, k);
    if(i>= (this->domain->get_num_grid_1())-1)
        point_data_r = new GridPoint(0., 0., std::complex<double>{0,0});
    auto point_data_u = this-> domain->at(i, j+1, k);
    if(j >= (this->domain->get_num_grid_2())-1)
        point_data_u = new GridPoint(0., 0., std::complex<double>{0,0});
    
        
    //potential at x, y 
    double V_ij = this->potential_func(point_data->x, point_data->y);
    //g * |psi(x,y)|^2 
    double additional_term = (this-> g) * (std::abs(point_data->wave_function))*(std::abs(point_data->wave_function));
    //Set infinitesimal value 
    double dx = this->domain->get_infinitesimal_distance1();
    double dy = this->domain->get_infinitesimal_distance2();
    //df denote time differential of dt (d(psi)/dt)
    // = (laplace - V-g|psi|^2) psi
    std::complex<double> df = 
        +((point_data_r->wave_function)+(point_data_l->wave_function)-(point_data->wave_function) * std::complex<double>{2})/ (std::complex<double>{dx*dx})
        +((point_data_u->wave_function)+(point_data_d->wave_function)-(point_data->wave_function) * std::complex<double>{2})/ (std::complex<double>{dy*dy})
        - (V_ij+additional_term) * (point_data->wave_function);
    df *= std::complex<double>{0,1}; 
    return df; 
};

/**
 * @brief Update k+1 th  wave function value with k th values 
 * 
 * @param k 
 */
void ForwardEulerRectangularSolver::solve_single_time(int k)
{   int Nx = this->domain->get_num_grid_1();
    int Ny = this->domain->get_num_grid_2();
    double dt = this->domain->get_dt();
    for(int i=0; i<Nx; ++i){
        for(int j=0; j<Ny; ++j){
            (this->domain->at(i, j, k+1)->wave_function) = this->domain->at(i, j, k)->wave_function + dt * this->temporal_equation(i,j,k);
        }
    }
}
void ForwardEulerRectangularSolver::applyInitialCondition(){
    this -> initialCondition.assign_to_domain(this->domain);
}
void ForwardEulerRectangularSolver::solve(){
    this->applyInitialCondition();
    int time_length = this->domain->get_num_times();
    for(int k=0; k<time_length-1; ++k){
        solve_single_time(k);
    }
}


