#pragma once 
#include <vector> 
#include <map>
#include <string.h>
using namespace std;

class BaseSweeper{
    public: 
        BaseSweeper() = default;
        BaseSweeper(float start, float end, int num, bool endpoint=false);
        float get_value_from_idx(int idx);
        float get_start(); 
        float get_end(); 
        int get_number_of_pts();
        //TODO
        //void set_info(map<string, string> info_map);
        virtual void get_MPI_info(int rank, int size);
        virtual void get_CUDA_info();
        ~BaseSweeper();

    protected: 
        float start;
        float end; 
        int num; 
        bool endpoint; 
        vector<float> num_list; 
        void generate_num_list();

        //MPI info 
        int rank=0 ;
        int size=0;
        bool MPI_use;
        bool CUDA_use;
        
};