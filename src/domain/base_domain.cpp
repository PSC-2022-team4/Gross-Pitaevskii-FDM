#include "base_domain.h"
#include <iostream>

GridPoint::GridPoint(double x_, double y_, std::complex<double> wave_function_) : x(x_), y(y_), wave_function(wave_function_){}

BaseSpatialGrid::BaseSpatialGrid(int num_grid_1, int num_grid_2)
{
    this->num_grid_1 = num_grid_1;
    this->num_grid_2 = num_grid_2;
    this->spatial_data = std::vector<std::vector<GridPoint>>(num_grid_1);
    for (auto i = 0; i < num_grid_1; ++i)
    {
        this->spatial_data[i] = std::vector<GridPoint>(num_grid_2);
    }
}
double BaseSpatialGrid::get_infinitesimal_distance1(){
    return this->infinitesimal_distance_1;
}
double BaseSpatialGrid::get_infinitesimal_distance2(){
    return this->infinitesimal_distance_2;
}
GridPoint * BaseSpatialGrid::at(int index_1, int index_2){
    return &this->spatial_data[index_1%this->num_grid_1][index_2%this->num_grid_2];
}
BaseDomain::BaseDomain(
    int num_grid_1,
    int num_grid_2,
    double t_start,
    double t_end,
    int num_times)
{
    this->t_start = t_start;
    this->t_end = t_end;
    this->num_times = num_times;
    this->dt = (t_end - t_start) / (num_times - 1);
    this->domain_data = std::vector<BaseSpatialGrid>(num_times);
    this->times = std::vector<double>(num_times);
    this->num_grid_1 = num_grid_1;
    this->num_grid_2 = num_grid_2;
    for (auto i = 0; i < num_times; ++i)
    {
        this->domain_data[i] = BaseSpatialGrid(num_grid_1, num_grid_2);
        this->times[i] = t_start + dt * i;
    }
}
double BaseDomain::get_t_start(){
    return this->t_start;
}

double BaseDomain::get_t_end(){
    return this->t_end;
}
double BaseDomain::get_dt(){
    return this->dt;
}
int BaseDomain::get_num_times()
{
    return this->num_times;
}

int BaseDomain::get_num_grid_1(){
    return this->num_grid_1;
}
int BaseDomain::get_num_grid_2(){
    return this->num_grid_2;
}

GridPoint * BaseDomain::at(int index_1, int index_2, int time_index){
    return this->domain_data[time_index].at(index_1, index_2);
}

void BaseDomain::assign_initial_value(int index_1, int index_2, std::complex<double> value){
    this->at(index_1, index_2, 0)->wave_function = value;
}
void BaseDomain::assign_wave_function(int index_1, int index_2, int time_index, std::complex<double> value){
    this->at(index_1, index_2, time_index)->wave_function = value;
}

double BaseDomain::time_at(int time_index){
    return this->times[time_index];
}

double BaseDomain::get_infinitesimal_distance1(){
    return this->domain_data[0].get_infinitesimal_distance1();
}
double BaseDomain::get_infinitesimal_distance2(){
    return this->domain_data[0].get_infinitesimal_distance2();
}
