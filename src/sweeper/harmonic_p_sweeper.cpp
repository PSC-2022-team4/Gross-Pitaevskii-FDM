/**
 * @file harmonic_p_sweeper.cpp
 * @author Minyoung Kim, Gyeonghun Kim
 * @brief Implementation of methods in the harmonic potential strength sweeper class.
 * @version 0.1
 * @date 2022-06-10
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include "./harmonic_p_sweeper.h"
#include "../solver/serial_solver/crank_nicolson/cn_rect_solver.h"

#include "../../src/solver/parallel_solver/crank_nicolson/cn_rect_psolver.cuh"
#include <mpi.h>
/**
 * @brief Construct a new HPSweeper::HPSweeper object
 * 
 * @param start 
 * @param end 
 * @param num 
 * @param endpoint 
 */
HPSweeper::HPSweeper(float start, float end, int num, bool endpoint)
: BaseSweeper(start, end, num, endpoint){
    //numlist contains unisotropy of angular frequency
    //But, we can rescale the value with domain length. -> set omega_x = 1  
    //Hence unisotropy =  omega_y / omega_x = omega_y 
}
/**
 * @brief solve the equation with various parameters 
 * 
 * @param domain domain to solve
 * @param initial_condition 
 * @param g 
 */
void HPSweeper::run(RectangularDomain *domain, InitialCondition *initial_condition, float g ){

    //If MPI_use = false, CUDA_use = false, solve the equation serially 
    if (!(this -> MPI_use) && !(this -> CUDA_use)){
        if(print_info){
            cout<< "Running serially started"<<endl;
        }
        
        for(int i=0; i<this -> num ; ++i){
            //Set conditions 
            initial_condition->assign_to_domain(domain);
            HarmonicPotential * potential = new HarmonicPotential(1, get_value_from_idx(i));
            potential->calcualte_potential_in_grid(domain);
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
    //If MPI_use = true, CUDA_use = false, solve the equation with various processors.
    else if ((this->MPI_use) && !(this->CUDA_use)){
        if(print_info){
            cout<< "Running with only MPI started"<<endl;
        }
        //If number of tasks exceed number of processors, abort. 
        if (num > size){
            MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE);
        }
        //If my rank is less than the task number, solve the equation
        if (rank < num){
            initial_condition->assign_to_domain(domain);
            HarmonicPotential * potential =new HarmonicPotential(1, get_value_from_idx(rank));
            potential -> calcualte_potential_in_grid(domain);
            CNRectSolver* solver =new CNRectSolver(g, domain);
            solver->solve(1e-11, 101, "MPI_"+to_string(rank), this -> print_info, this->save_data);
            delete solver;
            delete potential;
        }else{
            // No job for extra processors 
            ;
        }
    }//If MPI_use = false and CUDA_use =true, solve the equation with cuda, serially
    else if (!(this->MPI_use) && (this->CUDA_use)){

        if(print_info){
            cout<< "Running with CUDA serially started"<<endl;
        }
        for(int i=0; i<this -> num ; ++i){
            //Set conditions 
            initial_condition->assign_to_domain(domain);
            HarmonicPotential * potential = new HarmonicPotential(1, get_value_from_idx(i));
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
    //If MPI_use = true, CUDA_use = true, solve the equation with N parameters simultaneously. 
    //N =(number of device) * (maximum solver in one device)
    
    else{
        if(print_info){
            cout<< "Running with CUDA & MPI started"<<endl;
        }
       //If number of tasks exceed number of processors, abort. 
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
                CNRectPSolver* solver =new CNRectPSolver(g, domain, casted_GPU_num);
                solver->solve(1e-11, 101, "MPI&CUDA_"+to_string(rank), this -> print_info, this->save_data);
                delete solver;
                delete potential;
            }
        }


        
    }
}