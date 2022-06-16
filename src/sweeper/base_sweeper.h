/**
 * @file base_sweeper.h
 * @author Minyoung Kim, Gyeonghun Kim
 * @brief Header for Base sweeper class 
 * @version 0.1
 * @date 2022-06-10
 * 
 * @copyright Copyright (c) 2022
 * 
 */

  
#pragma once 
#include <vector> 
#include <map>
#include <string.h>
using namespace std;
/**
 * @brief Base sweeper class 
 * 
 */
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
        virtual void set_MPI_info(int rank, int size);
        virtual void set_CUDA_info(int gpu_num, int gpu_limit);
        void set_print_info(bool print_info);
        void set_save_data(bool save_data);
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
        //CUDA info 
        bool CUDA_use;
        int gpu_num=1;
        int gpu_limit=3;
        //Save results 
        bool print_info=true; 
        bool save_data=true;     
};