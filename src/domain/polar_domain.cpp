# include "polar_domain.h"
#include <cmath>

/**
 * @brief Construct a new Polar Spatial Grid:: Polar Spatial Grid object
 *
 * @param num_grid_1 number of grids in radial axis
 * @param num_grid_2 number of grids in angular axis
 * @param r_start_ start radius
 * @param r_end_ end radius
 */
PolarSpatialGrid::PolarSpatialGrid(int num_grid_1, int num_grid_2, double r_start_, double r_end_)
    : BaseSpatialGrid(num_grid_1, num_grid_2), r_start(r_start_), r_end(r_end_)
{   
    double r = r_start;
    double dr = (r_start - r_end) / (num_grid_1 - 1);
    double theta;
    double dtheta = (2 * M_PI) / (num_grid_2);
    for (auto i = 0; i < num_grid_1; ++i)
    {
        theta = 0;
        for (auto j = 0; j < num_grid_2; ++j)
        {
            this->spatial_data[i][j] = GridPoint(r * std::cos(theta), r * std::sin(theta), std::complex<double> {0, 0});
            theta += dtheta;
        }
        r += dr;
    }
}
/**
 * @brief Construct a new Polar Domain:: Polar Domain object
 *
 * @param num_grid_1 number of grids in radial axis
 * @param num_grid_2 number of grids in angular axis
 * @param t_start start time
 * @param t_end end time
 * @param num_times number of time grids
 * @param r_start start radius
 * @param r_end end radius
 */
PolarDomain::PolarDomain(int num_grid_1, int num_grid_2, double t_start, double t_end, int num_times, double r_start, double r_end)
 : BaseDomain(num_grid_1, num_grid_2, t_start, t_end, num_times)
 {
     for (auto i = 0; i < num_times; ++i)
     {
         this->domain_data[i] = PolarSpatialGrid(num_grid_1, num_grid_2, r_start, r_end);
     }
}
