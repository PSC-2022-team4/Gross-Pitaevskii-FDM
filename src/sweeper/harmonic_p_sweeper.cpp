#include "./harmonic_p_sweeper.h"
#include "../solver/serial_solver/crank_nicolson/cn_rect_solver.h"

#include "../../src/solver/parallel_solver/crank_nicolson/cn_rect_psolver.cuh"
#include <mpi.h>

HPSweeper::HPSweeper(float start, float end, int num, bool endpoint)
: BaseSweeper(start, end, num, endpoint){
    //numlist contains anharmonicity of angular frequency 
}

void HPSweeper::run(RectangularDomain *domain, InitialCondition *initial_condition, float g ){


    if (!(this -> MPI_use) && !(this -> CUDA_use)){
        if(print_info){
            cout<< "Running serially started"<<endl;
        }
        
        for(int i=0; i<this -> num ; ++i){
            //Set conditions 
            initial_condition->assign_to_domain(domain);
            auto potential = new HarmonicPotential(1, get_value_from_idx(i));
            potential->calcualte_potential_in_grid(domain);
            g = this  -> get_value_from_idx(i);
            //Apply on solver
            CNRectSolver* solver =new CNRectSolver(g, domain);
            //solve using solver. It automatically save data
            solver->solve(1e-11, 101, to_string(i) , this-> print_info, this -> save_data);
            //reset domain to use in next iteration 
            domain->reset();
            delete solver;
            delete potential;
        }
    }
    else if ((this->MPI_use) && !(this->CUDA_use)){
        //If number of tasks exceed number of processors, abort. 
        if(print_info){
            cout<< "Running with only MPI started"<<endl;
        }
        if (num > size){
            MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE);
        }
        float g = 0 ; 
        if (rank < num){
            initial_condition->assign_to_domain(domain);
            auto potential =new HarmonicPotential(1, get_value_from_idx(rank));
            potential->calcualte_potential_in_grid(domain);
            CNRectSolver* solver =new CNRectSolver(g, domain);
            solver->solve(1e-11, 101, "MPI_"+to_string(rank), this -> print_info, this->save_data);
            delete solver;
            delete potential;
        }else{
            // No job for extra processors 
            ;
        }
    }else if (!(this->MPI_use) && (this->CUDA_use)){

        if(print_info){
            cout<< "Running with CUDA serially started"<<endl;
        }
        for(int i=0; i<this -> num ; ++i){
            //Set conditions 
            initial_condition->assign_to_domain(domain);
            auto potential = new HarmonicPotential(1, get_value_from_idx(i));
            potential->calcualte_potential_in_grid(domain);
            //Apply on solver
            CNRectPSolver* solver =new CNRectPSolver(g, domain, 0);//solve using solver. It automatically save data
            solver->solve(1e-11, 101, "CUDA_"+to_string(i) , this-> print_info, this -> save_data);
            //reset domain to use in next iteration 
            domain->reset();
            delete solver;
            delete potential;
        }
        
    }
    else{
        if(print_info){
            cout<< "Running with CUDA & MPI started"<<endl;
        }
        float g=0;
        
        if (num > size){
            MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE);
        }
        int max_pallel_tasks = this -> gpu_num * this ->gpu_limit;
        int repeat = num / max_pallel_tasks+(num%max_pallel_tasks !=0);
        for(int i=0; i<repeat; ++i){
            if (i*repeat < max_pallel_tasks && rank < (i+1)*max_pallel_tasks && rank < num ){
                int casted_GPU_num = (rank % max_pallel_tasks) % gpu_limit;
                initial_condition->assign_to_domain(domain);
                auto potential = new HarmonicPotential(1, get_value_from_idx(rank));
                potential->calcualte_potential_in_grid(domain);
                //g = this  -> get_value_from_idx(rank);
                CNRectPSolver* solver =new CNRectPSolver(g, domain, casted_GPU_num);
                solver->solve(1e-11, 101, "MPI&CUDA_"+to_string(rank), this -> print_info, this->save_data);
                delete solver;
                delete potential;
            }
        }


        
    }
}