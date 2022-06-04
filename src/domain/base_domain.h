#pragma once
#include <complex>
#include <vector>

class GridPoint
{
    public:
        GridPoint();
        GridPoint(double x, double y, std::complex<double> wave_function);
        double x, y;
        std::complex<double> wave_function;
};

class BaseSpatialGrid
{
    public:
        BaseSpatialGrid();
        BaseSpatialGrid(int num_grid_1, int num_grid_2);
        GridPoint * at(int index_1, int index_2);
        double get_infinitesimal_distance1();
        double get_infinitesimal_distance2();
        

    protected:
        std::vector<std::vector<GridPoint>> spatial_data;
        double infinitesimal_distance_1, infinitesimal_distance_2;
        int num_grid_1, num_grid_2;
};

class BaseDomain
{
public:
    BaseDomain();
    BaseDomain(int num_grid_1, int num_grid_2, double t_start, double t_end, int num_times);
    double get_t_start();
    double get_t_end();
    double get_dt();
    int get_num_times();
    int get_num_grid_1();
    int get_num_grid_2();
    GridPoint * at(int index_1, int index_2, int time_index);
    void assign_initial_value(int index_1, int index_2, std::complex<double> value);
    double time_at(int time_index);
    double get_infinitesimal_distance1();
    double get_infinitesimal_distance2();

protected:
    std::vector<BaseSpatialGrid> domain_data;
    std::vector<double> times;
    double t_start, t_end, dt;
    int num_times, num_grid_1, num_grid_2;
};