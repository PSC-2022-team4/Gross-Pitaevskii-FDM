#include "base_domain.h"
#include <iostream>

GridPoint::GridPoint(double x_, double y_, std::complex<double> wave_function_) : x(x_), y(y_), wave_function(wave_function_){}
// BaseSpatialGrid::BaseSpatialGrid(int num_grid_1, int num_grid_2)
// {
//     this->num_grid_1 = num_grid_1;
//     this->num_grid_2 = num_grid_2;
//     this->spatial_data.resize(num_grid_1);
//     for (auto i = 0; i < num_grid_1; ++i)
//     {
//         this->spatial_data[i].resize(num_grid_2);
//     }
// }

// GridPoint BaseSpatialGrid::at(int index_1, int index_2){
//     return this->spatial_data[index_1][index_2];
// }

// BaseDomain::BaseDomain(
//     int num_grid_1, 
//     int num_grid_2,
//     double t_start,
//     double t_end, 
//     int num_times)
// {
//     this->t_start = t_start;
//     this->t_end = t_end;
//     this->num_times = num_times;
//     this->dt = (t_end - t_start) / (num_times - 1);
//     this->domain_data = std::vector<BaseSpatialGrid>(num_times);
//     std::cout << num_times << std::endl;
//     for (auto i = 0; i < num_times; ++i)
//     {
//         this->domain_data[i] = BaseSpatialGrid(num_grid_1, num_grid_2);
//         this->times[i] = t_start + dt * i;
//     }
// }
// double BaseDomain::get_t_start(){
//     return this->t_start;
// }
// double BaseDomain::get_t_end(){

//     return this->t_end;
// }
// double BaseDomain::get_dt(){

//     return this->dt;
// }
// int BaseDomain::get_num_times()
// {
//     return this->num_times;
// }