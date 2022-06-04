#include "base_domain.h"
#include <iostream>

<<<<<<< HEAD
GridPoint::GridPoint(double x_, double y_, std::complex<double> wave_function_) : x(x_), y(y_), wave_function(wave_function_){}

=======
GridPoint::GridPoint(){}
/**
 * @brief Construct a new Grid Point:: Grid Point object
 * 
 * @param x_ grid point x  
 * @param y_ grid point y 
 * @param wave_function_ wave function value at x, y 
 */
GridPoint::GridPoint(double x_, double y_, std::complex<double> wave_function_) : x(x_), y(y_), wave_function(wave_function_){}

BaseSpatialGrid::BaseSpatialGrid(){}
/**
 * @brief Construct a new Base Spatial Grid:: Base Spatial Grid object
 * 
 * @param num_grid_1 Number of grid of first axis. index: 0 to num_grid_1-1
 * @param num_grid_2 Number of grid of second axis. index: 0 to num_grid_2-1
 */
>>>>>>> ae9304929e6d3ef7fecba2d21225c249a7c3bd33
BaseSpatialGrid::BaseSpatialGrid(int num_grid_1, int num_grid_2)
{
    this->num_grid_1 = num_grid_1;
    this->num_grid_2 = num_grid_2;
    //Spatial_data: 2D matrix of grid points 
    this->spatial_data = std::vector<std::vector<GridPoint>>(num_grid_1);
    for (auto i = 0; i < num_grid_1; ++i)
    {
        this->spatial_data[i] = std::vector<GridPoint>(num_grid_2);
    }
}
//Getter function of infinitesimal_distance1
double BaseSpatialGrid::get_infinitesimal_distance1(){
    return this->infinitesimal_distance_1;
}

//Getter function of infinitesimal_distance2
double BaseSpatialGrid::get_infinitesimal_distance2(){
    return this->infinitesimal_distance_2;
}
//Getter function of grid point at index_1, index_2.
//roll up for boundary  
GridPoint * BaseSpatialGrid::at(int index_1, int index_2){
    return &this->spatial_data[index_1%this->num_grid_1][index_2%this->num_grid_2];
}
<<<<<<< HEAD
=======

BaseDomain::BaseDomain(){}
/**
 * @brief Construct a new Base Domain:: Base Domain object
 *        It contains spatial and temporal domain 
 * 
 * @param num_grid_1 Number of grid of first axis. index: 0 to num_grid_1-1
 * @param num_grid_2 Number of grid of second axis. index: 0 to num_grid_2-1
 * @param t_start start time 
 * @param t_end end time 
 * @param num_times number of time points. 
 */
>>>>>>> ae9304929e6d3ef7fecba2d21225c249a7c3bd33
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
    //dt: interval of time 
    this->dt = (t_end - t_start) / (num_times - 1);
    // vector of spatial domain. domain_data[i] : spatial domain at time i * dt
    this->domain_data = std::vector<BaseSpatialGrid>(num_times);
    //vector of time list 
    this->times = std::vector<double>(num_times);
    this->num_grid_1 = num_grid_1;
    this->num_grid_2 = num_grid_2;
    for (auto i = 0; i < num_times; ++i)
    {
        this->domain_data[i] = BaseSpatialGrid(num_grid_1, num_grid_2);
        this->times[i] = t_start + dt * i;
    }
}
//Getter functions 
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

double BaseDomain::get_infinitesimal_distance1(){
    return this->domain_data[0].get_infinitesimal_distance1();
}
double BaseDomain::get_infinitesimal_distance2(){
    return this->domain_data[0].get_infinitesimal_distance2();
}
//Getter function of grid point at index_1, index_2 and time_index 
GridPoint * BaseDomain::at(int index_1, int index_2, int time_index){
    return this->domain_data[time_index].at(index_1, index_2);
}

/**
 * @brief Assign initial value 
 * 
 * @param index_1 x = x_start + index_1*dx
 * @param index_2 y = y_start + index_2 * dy 
 * @param value initial value at x, y 
 */
void BaseDomain::assign_initial_value(int index_1, int index_2, std::complex<double> value){
    this->at(index_1, index_2, 0)->wave_function = value;
}
void BaseDomain::assign_wave_function(int index_1, int index_2, int time_index, std::complex<double> value){
    this->at(index_1, index_2, time_index)->wave_function = value;
}
//Get time using time index 
double BaseDomain::time_at(int time_index){
    return this->times[time_index];
}
