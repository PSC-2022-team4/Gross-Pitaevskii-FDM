#include "src/domain/polar_domain.cpp"
#include "src/utils.h"
#include <iostream>


bool test_polar_spatial_grid()
{
    auto spatial_grid = PolarSpatialGrid(10, 10, 1, 2);
    bool all_passed = true;
    if (!is_close(spatial_grid.at(0, 0).x, 1., 1e-12))
    {
        all_passed = false;
    }
    if (!is_close(spatial_grid.at(0, 0).y, 0., 1e-12))
    {
        all_passed = false;
    }
    if (!is_close(spatial_grid.at(10, 10).wave_function.real(), 0., 1e-12))
    {
        all_passed = false;
    }
    if (!is_close(spatial_grid.at(10, 10).wave_function.imag(), 0., 1e-12))
    {
        all_passed = false;
    }

    return all_passed;
}

bool test_polar_domain_contructor()
{
    auto domain = PolarDomain(10, 10, 0, 10, 11, 1, 2);
    bool all_passed = true;
    if (!is_close(domain.get_t_start(), 0., 1e-12))
    {
        all_passed = false;
    }
    if (!is_close(domain.get_t_end(), 10., 1e-12))
    {
        all_passed = false;
    }
    if (domain.get_num_times() != 11)
    {
        all_passed = false;
    }
    if (!is_close(domain.get_dt(), 1., 1e-12))
    {
        all_passed = false;
    }
    return all_passed;
}

bool test_all_polar_domain()
{
    if (test_polar_spatial_grid())
    {
        std::cout << "test_polar_spatial_grid succeeded!" << std::endl;
    }
    else
    {
        std::cout << "test_polar_spatial_grid constuctor failed!" << std::endl;
    }
    if (test_polar_domain_contructor())
    {
        std::cout << "test_polar_domain_contructor succeeded!" << std::endl;
    }
    else
    {
        std::cout << "test_polar_domain_contructor failed!" << std::endl;
    }
}