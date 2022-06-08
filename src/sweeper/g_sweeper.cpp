#include "g_sweeper.h"
#include "../solver/serial_solver/crank_nicolson/cn_rect_solver.h"

gSweeper::gSweeper(float start, float end, int num, bool endpoint, bool MPI_use, bool CUDA_use)
: BaseSweeper(start, end, num, endpoint){
    
}

void gSweeper::run(RectangularDomain *domain,InitialCondition *initial_condition, HarmonicPotential *potential ){


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
            //solve using solver it 
            solver->solve(1e-11, 101, "test_"+to_string(i));
            domain->reset();
            delete solver;
        }

    }



}