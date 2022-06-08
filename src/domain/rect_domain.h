#pragma once
#include "base_domain.h"

class RectangularSpatialGrid : public BaseSpatialGrid
{
public:
    RectangularSpatialGrid() = default;
    RectangularSpatialGrid(int num_grid_1, int num_grid_2, float x_start, float x_end, float y_start, float y_end);
    ~RectangularSpatialGrid();

private:
    float x_start;
    float x_end;
    float y_start;
    float y_end;
};

class RectangularDomain : public BaseDomain
{
public:
    RectangularDomain() = default;
    RectangularDomain(int num_grid_1, int num_grid_2, float t_start, float t_end, int num_times, float x_start, float x_end, float y_start, float y_end);
    float get_x_start();
    float get_y_start();
    float get_x_end();
    float get_y_end();
    void generate_single_txt_file(std::string filename);
    ~RectangularDomain();
    // a grid to save potential values
    BaseSpatialGrid *potential_grid;
    void update_time();

private:
    float x_start;
    float x_end;
    float y_start;
    float y_end;
};