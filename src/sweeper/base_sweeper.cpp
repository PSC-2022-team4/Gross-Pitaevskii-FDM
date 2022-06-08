#include "base_sweeper.h"
#include <cassert>
using namespace std;

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

void BaseSweeper::get_MPI_info(int rank, int size){
    this -> rank = rank; 
    this -> size = size;
    this -> MPI_use=true;
}

void BaseSweeper::get_CUDA_info(){
    this -> CUDA_use=true;
}
BaseSweeper::~BaseSweeper(){};

//Set functions 

// void BaseSweeper::set_info(std::map<std::string, std::string> info_map){
//     //string solver_type = info_map["solver"];
//     //TODO 

// };
