/**
 * @file base_sweeper.cpp
 * @author Minyoung Kim, Gyeonghun Kim
 * @brief Implementation file for base sweeper
 * @version 0.1
 * @date 2022-06-05
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "base_sweeper.h"
#include <cassert>
using namespace std;
/**
 * @brief Construct a new Base Sweeper:: Base Sweeper object
 * 
 * @param start_ start value of paramters
 * @param end_ end value of parameters
 * @param num_ number of parameters
 * @param endpoint_ bool type. If true, sweep until end_, else, sweep until just before end_
 */
BaseSweeper::BaseSweeper(float start_, float end_, int num_, bool endpoint_):
start(start_), end(end_), num(num_), endpoint(endpoint_){
    
    this->generate_num_list(); 

}
/**
 * @brief generate number of points that sweep from start to end 
 * 
 */
void BaseSweeper::generate_num_list(){
    this->num_list = std::vector<float>(this->num);
    float d = 0.;
    if(this->endpoint){
        d = (this -> end - this-> start ) / float(this-> num -1 );
    }else{
        d = (this -> end - this-> start ) / float(this-> num);
    }

    for (int i=0; i<this -> num; ++i){
        this->num_list[i] = this ->start + d*i; 
    }

}

//Getter functions 
float BaseSweeper::get_value_from_idx(int idx){
    return this->num_list[idx];
}

float BaseSweeper::get_start(){
    return start; 
}
float BaseSweeper::get_end(){
    return end; 
}
int BaseSweeper::get_number_of_pts(){
    return num; 
}
/**
 * @brief Setter function with MPI
 *        If information is update, set MPI_use to true to indicate MPI usage. 
 * 
 * @param rank current rank of processor
 * @param size total number of processors
 */
void BaseSweeper::set_MPI_info(int rank, int size){
    this -> rank = rank; 
    this -> size = size;
    this -> MPI_use=true;
}
/**
 * @brief Setter function with CUDA
 *        If information is update, set CUDA_use to true to indicate CUDA usage.
 * 
 * @param gpu_num total number of GPU devices
 * @param gpu_limit maximum number of solvers using one GPU device 
 */
void BaseSweeper::set_CUDA_info(int gpu_num, int gpu_limit){
    this -> CUDA_use=true;
    this -> gpu_num = gpu_num;
    this -> gpu_limit = gpu_limit;
}
/**
 * @brief Setter function of print_info
 * 
 * @param print_info 
 */
void BaseSweeper::set_print_info(bool print_info){
    this -> print_info = print_info;
}
/**
 * @brief Setter function of save_data
 * 
 * @param save_data 
 */
void BaseSweeper::set_save_data(bool save_data){
    this -> save_data = save_data;
}
/**
 * @brief Destroy the Base Sweeper:: Base Sweeper object
 * 
 */
BaseSweeper::~BaseSweeper(){};

