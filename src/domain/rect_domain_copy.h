#pragma once
#include "base_domain_copy.h"


class RectangularSpatialGrid: public BaseSpatialGrid
{
public:
    RectangularSpatialGrid() = default;
    RectangularSpatialGrid(int num_grid_1, int num_grid_2, double x_start, double x_end, double y_start, double y_end);

private: 
    double x_start; 
    double x_end; 
    double y_start; 
    double y_end;
};

class RectangularDomain : public BaseDomain
{
public:
    RectangularDomain() = default;
    RectangularDomain(int num_grid_1, int num_grid_2, double t_start, double t_end, int num_times, double x_start, double x_end, double y_start, double y_end);
    double get_x_start();
    double get_y_start();
    double get_x_end();
    double get_y_end();
private:
    double x_start; 
    double x_end; 
    double y_start; 
    double y_end;
};