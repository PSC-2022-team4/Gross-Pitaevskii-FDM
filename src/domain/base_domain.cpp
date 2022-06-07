#include "base_domain.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <sstream>
#include <experimental/filesystem>

namespace fs = std::experimental::filesystem;

GridPoint::GridPoint(double x_, double y_, std::complex<double> wave_function_) : x(x_), y(y_), value(wave_function_){}


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

void BaseSpatialGrid::normalize(){
    double sum = 0.;
    for (auto i = 0; i < this->num_grid_1; ++i){
        for (auto j = 0; j < this->num_grid_2;++j){
            auto wave_func = this->at(i, j)->value;
            sum += std::pow(std::abs(wave_func), 2);
        }
    }
    sum = std::sqrt(sum * this->infinitesimal_distance_1 * this->infinitesimal_distance_2);
    for (auto i = 0; i < this->num_grid_1; ++i)
    {
        for (auto j = 0; j < this->num_grid_2; ++j)
        {
            this->at(i, j)->value /= sum;
        }
    }
}
// Getter function of infinitesimal_distance1
double BaseSpatialGrid::get_infinitesimal_distance1()
{
    return this->infinitesimal_distance_1;
}

// Getter function of infinitesimal_distance2
double BaseSpatialGrid::get_infinitesimal_distance2()
{
    return this->infinitesimal_distance_2;
}
// Getter function of grid point at index_1, index_2.
// roll up for boundary
GridPoint *BaseSpatialGrid::at(int index_1, int index_2)
{
    return &this->spatial_data[index_1 % this->num_grid_1][index_2 % this->num_grid_2];
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
    // dt: interval of time
    this->dt = (t_end - t_start) / (num_times - 1);
    // vector of spatial domain. domain_data[i] : spatial domain at time i * dt
    this->domain_data = std::vector<BaseSpatialGrid>(num_times);
    // vector of time list
    this->times = std::vector<double>(num_times);
    this->num_grid_1 = num_grid_1;
    this->num_grid_2 = num_grid_2;
    for (auto i = 0; i < num_times; ++i)
    {
        this->domain_data[i] = BaseSpatialGrid(num_grid_1, num_grid_2);
        this->times[i] = t_start + dt * i;
    }
}

void BaseDomain::normalize(int time_index){
    this->domain_data[time_index].normalize();
}
// Getter functions
double BaseDomain::get_t_start()
{
    return this->t_start;
}

double BaseDomain::get_t_end()
{
    return this->t_end;
}
double BaseDomain::get_dt()
{
    return this->dt;
}
int BaseDomain::get_num_times()
{
    return this->num_times;
}

int BaseDomain::get_num_grid_1()
{
    return this->num_grid_1;
}
int BaseDomain::get_num_grid_2()
{
    return this->num_grid_2;
}

double BaseDomain::get_infinitesimal_distance1()
{
    return this->domain_data[0].get_infinitesimal_distance1();
}
double BaseDomain::get_infinitesimal_distance2()
{
    return this->domain_data[0].get_infinitesimal_distance2();
}
// Getter function of grid point at index_1, index_2 and time_index
GridPoint *BaseDomain::at(int index_1, int index_2, int time_index)
{
    return this->domain_data[time_index].at(index_1, index_2);
}

/**
 * @brief Assign initial value
 *
 * @param index_1 x = x_start + index_1*dx
 * @param index_2 y = y_start + index_2 * dy
 * @param value initial value at x, y
 */
void BaseDomain::assign_initial_value(int index_1, int index_2, std::complex<double> value)
{
    this->at(index_1, index_2, 0)->value = value;
}
void BaseDomain::assign_wave_function(int index_1, int index_2, int time_index, std::complex<double> value)
{
    this->at(index_1, index_2, time_index)->value = value;
}
// Get time using time index
double BaseDomain::time_at(int time_index)
{
    return this->times[time_index];
}
/**
 * @brief create directory and return directory name with "/"
 *         directory name : "./results/%d-%m-%Y %H-%M-%S_"+info
 *         for exmaple, "./results/"
 * @param info  information about domain type and solver type
 * @return std::string directory name with "/"
 */
std::string BaseDomain::generate_directory_name(std::string info)
{
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);

    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%d-%H-%M-%S");
    auto str = oss.str();
    std::cout << str << std::endl;
    std::string directory_name = "../results/" + str + "_" + info;

    if (fs::create_directory(directory_name.c_str()))
    {
        std::cout << "Created directory " << directory_name << std::endl;
    }
    else
    {
        std::cout << "Creating directory failed" << std::endl;
    };
    return directory_name + "/";
}

/**
 * @brief export results of probability(|psi|^2) as txt
 *        it generates (num_times) txt files
 *        each txt file contains  |psi|^2 data
 * @param info To generate directory
 */
std::string BaseDomain::generate_txt_file(std::string info)
{
    std::string base_filename = this->generate_directory_name(info);
    std::string filename = "";
    for (int t = 0; t < this->num_times; ++t)
    {
        filename = std::string("probability_") + std::to_string(t);
        filename = base_filename + filename;
        generate_single_txt_file(&domain_data[t], filename);
    }
    std::cout << this->num_times;
    std::cout << " text files are generated in \n";
    std::cout << base_filename << std::endl;

    return base_filename;
}
/**
 * @brief write txt file at certain t
 *
 * @param grid grid at t
 * @param filename file name to store data
 */
void BaseDomain::generate_single_txt_file(BaseSpatialGrid * grid, std::string filename)
{
    std::ofstream outfile(filename + ".txt");
    outfile << "x, y, probability" << std::endl;
    for (auto i = 0; i < num_grid_1; ++i)
    {
        for (auto j = 0; j < num_grid_2; ++j)
        {
            double magnitude = std::abs(grid->at(i, j)->value);
            outfile << grid->at(i, j)->x << ", " << grid->at(i, j)->y << ", ";
            outfile << magnitude * magnitude;
            outfile << std::endl;
        }
    }
    outfile.close();
};
