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
        std::complex<float> wave_function;
        ~GridPoint();
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
        ~BaseSpatialGrid();

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
    //For boundary 
    GridPoint * get_null_gridpt();
    GridPoint * at(int index_1, int index_2, int time_index);
    void assign_initial_value(int index_1, int index_2, std::complex<float> value);
    void assign_wave_function(int index_1, int index_2, int time_index, std::complex<float> value);
    float time_at(int time_index);
    float get_infinitesimal_distance1();
    float get_infinitesimal_distance2();
    void generate_directory_name(std::string info);
    void generate_single_txt_file(std::string filename);

    void normalize(int time_index);
    void print_directory_info();
    ~BaseDomain();

protected:
    BaseSpatialGrid* old_grid; 
    BaseSpatialGrid* current_grid; 
    GridPoint * null_gridpt;
    float t_start, t_end, dt;
    int num_times, num_grid_1, num_grid_2;
    int current_time_index = 0;
    std::string PATH;
private:
    void update_time();
};