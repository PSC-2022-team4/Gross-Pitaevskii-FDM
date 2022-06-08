#pragma once
#include <complex>
#include <vector>
#include <string.h>
class GridPoint
{
    public:
        GridPoint() = default;
        GridPoint(float x, float y, std::complex<float> wave_function);
        float x, y;
        std::complex<float> value;
};

class BaseSpatialGrid
{
    public:
        BaseSpatialGrid() = default;
        BaseSpatialGrid(int num_grid_1, int num_grid_2);
        GridPoint * at(int index_1, int index_2);
        float get_infinitesimal_distance1();
        float get_infinitesimal_distance2();
        void normalize();

    protected:
        std::vector<std::vector<GridPoint>> spatial_data;
        float infinitesimal_distance_1, infinitesimal_distance_2;
        int num_grid_1, num_grid_2;
};

class BaseDomain
{
public:
    BaseDomain();
    BaseDomain(int num_grid_1, int num_grid_2, float t_start, float t_end, int num_times);
    float get_t_start();
    float get_t_end();
    float get_dt();
    int get_num_times();
    int get_num_grid_1();
    int get_num_grid_2();
    GridPoint * at(int index_1, int index_2, int time_index);
    void assign_initial_value(int index_1, int index_2, std::complex<float> value);
    void assign_wave_function(int index_1, int index_2, int time_index, std::complex<float> value);
    float time_at(int time_index);
    float get_infinitesimal_distance1();
    float get_infinitesimal_distance2();
    std::string generate_txt_file(std::string info);
    void normalize(int time_index);

protected:
    std::vector<BaseSpatialGrid> domain_data;
    std::vector<float> times;
    float t_start, t_end, dt;
    int num_times, num_grid_1, num_grid_2;
    std::string generate_directory_name(std::string info);
    void generate_single_txt_file(BaseSpatialGrid* grid, std::string filename);
};