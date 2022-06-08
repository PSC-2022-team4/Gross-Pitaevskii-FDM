#pragma once 
#include <vector> 

class BaseSweeper{
    public: 
        BaseSweeper() = default;
        BaseSweeper(float start, float end, int num, bool endpoint=false);
        float get_value_from_idx(int idx);
        float get_start(); 
        float get_end(); 
        int get_number_of_pts();
        ~BaseSweeper();

    private: 
        float start;
        float end; 
        int num; 
        bool endpoint; 
        std::vector<float> num_list; 
        void generate_num_list();


};