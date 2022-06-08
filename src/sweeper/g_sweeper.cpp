#include "./g_sweeper.h"
#include "../solver/serial_solver/crank_nicolson/cn_rect_solver.h"
#include <mpi.h>
gSweeper::gSweeper(float start, float end, int num, bool endpoint, bool MPI_use, bool CUDA_use)
: BaseSweeper(start, end, num, endpoint){
    
}

void gSweeper::run(RectangularDomain *domain, InitialCondition *initial_condition, HarmonicPotential *potential ){


    if (!(this -> MPI_use) && !(this -> CUDA_use)){
        cout<< "Running serially started"<<endl;
        float g=0;
        for(int i=0; i<this -> num ; ++i){
            //Set conditions 
            initial_condition->assign_to_domain(domain);
            potential->calcualte_potential_in_grid(domain);
            g = this  -> get_value_from_idx(i);
            //Apply on solver
            CNRectSolver* solver =new CNRectSolver(g, domain);
            //solve using solver. It automatically save data
            solver->solve(1e-11, 101, "_"+to_string(i));
            //reset domain to use in next iteration 
            domain->reset();
            delete solver;
        }
    }
    else if ((this->MPI_use) && !(this->CUDA_use)){
        //If number of tasks exceed number of processors, abort. 
        if (num > size){
            MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE);
        }
        float g = 0 ; 
        if (rank < num){
            initial_condition->assign_to_domain(domain);
            potential->calcualte_potential_in_grid(domain);
            g = this  -> get_value_from_idx(rank);
            CNRectSolver* solver =new CNRectSolver(g, domain);
            solver->solve(1e-11, 101, "_"+to_string(rank));
        }else{
            // No job for extra processors 
            ;
        }
    }else if (!(this->MPI_use) && (this->CUDA_use)){
        //Using only CUDA
        ;
    }
    else{
        //both MPI and CUDA 
        ;
    }
}