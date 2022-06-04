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
    // uto infinitesimal_distance_1 = this->domain.get
    auto point_data = this->domain->at(i, j, k);
    auto point_data_left = this->domain->at(i - 1, j, k);
    auto point_data_right = this->domain->at(i+1, j, k);
    auto point_data_up = this->domain->at(i, j-1, k);
    auto point_data_down = this->domain->at(i, j+1, k);




    // double this->potential()
    // std::vector<BaseSpatialGrid> * domain_data = &(this-> domain->domain_data) ;
};
