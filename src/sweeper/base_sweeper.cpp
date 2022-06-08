#include <base_sweeper.h>
#include <cassert>

BaseSweeper::BaseSweeoer(float start_, float end_, int num_, bool endpoint_):
start(start_), end(end_), num(num_), endpoint(endpoint_){
    
    this->generate_num_list(); 
}
/**
 * @brief generate number of points that sweep from start to end 
 * 
 */
void BaseSweeper::generate_num_list(){
    this->num_list = std::vector<float>(this->num);
    
    if(this->endpoint){
        float d = (this -> end - this-> start ) / float(this-> num -1 );
    }else{
        float d = (this -> end - this-> start ) / float(this-> num);
    }

    for (int i=0; i<this -> num; ++i){
        this->num_list[i] = start_ + d; 
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

BaseSweeper::~BaseSweeper(){};