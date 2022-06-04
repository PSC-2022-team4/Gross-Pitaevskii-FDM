#include "base_serial_solver.h"

/**
 * @brief Construct a new Base Solver:: Base Solver object
 *        Default constructor
 * 
 */
BaseSolver::BaseSolver(){};

BaseSolver::BaseSolver(InitialCondition initialCondition){
    
    this -> initialCondition = initialCondition; 

};
