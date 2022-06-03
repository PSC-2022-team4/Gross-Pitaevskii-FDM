#include "base_domain.h"

class RectangularSpatialGrid: private BaseSpatialGrid
{
public:
    RectangularSpatialGrid(int num_grid_1, int num_grid_2, double x_start, double x_end, double y_start, double y_end);
};

class RectangularDomain : private BaseDomain
{
public:
    RectangularDomain(int num_grid_1, int num_grid_2, double t_start, double t_end, int num_times);
};