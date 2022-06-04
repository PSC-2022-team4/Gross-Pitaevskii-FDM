#include "src/serial_solver/forward_euler/forward_euler_rectangular_solver.h"

ForwardEulerRectangularSolver::ForwardEulerRectangularSolver(
    InitialCondition initialCondition,
    std::function<double(double, double)> potential, 
    double g,
    RectangularDomain rectangularDomain)
    :BaseSolver(initialCondition, potential, g, rectangularDomain)
    {
        //Cast basedomain to rectangularDomain
        (RectangularDomain) *(this->domain);
    };

void ForwardEulerRectangularSolver::applyInitialCondition(){
    this -> initialCondition.assign_to_domain(this->domain);
}

/**
 * @brief Time differential of phi 
 * 
 * @param i index for x 
 * @param j index for y 
 * @param k index for time 
 * @return std::complex<double> 
 */
std::complex<double> ForwardEulerRectangularSolver::temporal_equation(int i, int j, int k){
    this->domain->at(i, j, k);   
    auto point_data = this-> domain->at(i, j, k);
    //double this->potential()
    //std::vector<BaseSpatialGrid> * domain_data = &(this-> domain->domain_data) ;

};


