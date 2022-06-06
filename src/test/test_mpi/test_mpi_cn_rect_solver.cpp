//#include "src/serial_solver/crank_nicolson/cn_rect_solver.h"

#include "src/serial_solver/crank_nicolson/cn_rect_solver_copy.h"
#include "src/utils.h"
#include <functional>
#include <iostream>
#include <complex>

bool test_mpi_cn_solver_creation(int rank, int size)
{
    bool all_passed = true;

    RectangularDomain *domain = (new RectangularDomain(41, 41, 0, 0.03, 3, -5, 5, -5, 5));
    /*
    Worked values :
    1) (41, 41, 0, 0.1, 101, -5, 5, -5, 5), initial_cond /9, potential 0.5, g=1
    2) (41, 41, 0, 1, 101, -5, 5, -5, 5), initial_cond /9, potential 0.5, g=1 only condensate 
    3) (41, 41, 0, 1, 101, -5, 5, -5, 5), initial_cond /9, potential 1, g=1
    4) (41, 41, 0, 1, 101, -5, 5, -5, 5), initial_cond 0(uniform), potential 1, g=-10 : start from uniform, it condensates 
    5) 
    */
    auto initial_cond_function = [](double x, double y)
    { //return std::complex<double>{1 * std::exp(-(x * x + y * y) *0/ (9))}; };
        return std::complex<double> {1.,0.};};
    auto *initial_condition = new InitialCondition(initial_cond_function);
    initial_condition->assign_to_domain(domain);

    auto potential = [rank](double x, double y)
    {
        return (double) 1 * (x * x + y * y);
    };
    double g = -10.;
    CNRectSolver solver = CNRectSolver(potential, g, domain);

    solver.solve(1e-13, 101, std::to_string(rank));

    return all_passed;
}

bool test_all_mpi_cn_solver(int rank, int size)
{
    if (test_mpi_cn_solver_creation(rank, size ))
    {
        std::cout << "test_crank_nicolson_solver_creation succeeded!" << std::endl;
    }
    else
    {
        std::cout << "test_crank_nicolson_solver_creation failed!" << std::endl;
    }
}
