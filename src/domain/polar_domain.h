#pragma once
#include "base_domain.h"
class PolarSpatialGrid : public BaseSpatialGrid
{
public:
    PolarSpatialGrid() = default;
    PolarSpatialGrid(int num_grid_1, int num_grid_2, double r_start_, double r_end_);

private:
    double r_start, r_end;
};

class PolarDomain : public BaseDomain
{
public:
    PolarDomain() = default;
    PolarDomain(int num_grid_1, int num_grid_2, double t_start, double t_end, int num_times, double r_start, double r_end);
};