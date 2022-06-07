#include "src/domain/rect_domain.h"
#include "src/utils.h"
#include <iostream>

bool test_rectangular_spatial_grid(){
    auto spatial_grid = RectangularSpatialGrid (101, 101, 0, 10, 0, 20);
    bool all_passed = true;
    if (!is_close(spatial_grid.at(10, 10)->x, 1., 1e-12))
    {
        all_passed = false;
    }
    if (!is_close(spatial_grid.at(10, 10)->y, 2., 1e-12))
    {
        all_passed = false;
    }
    if (!is_close(spatial_grid.at(10, 10)->value.real(), 0., 1e-12))
    {
        all_passed = false;
    }
    if (!is_close(spatial_grid.at(10, 10)->value.imag(), 0., 1e-12))
    {
        all_passed = false;
    }

    return all_passed;
}

bool test_rect_domain_contructor(){
    auto domain = RectangularDomain(100, 100, 0., 10., 11, 0, 10 , 0, 20 );
    bool all_passed = true;
    if (!is_close(domain.get_t_start(), 0., 1e-12)){
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

bool test_all_rect_domain()
{
    if (test_rectangular_spatial_grid())
    {
        std::cout << "Test rectangular spatial grid succeeded!" << std::endl;
    }
    else
    {
        std::cout << "Test rectangular spatial grid constuctor failed!" << std::endl;
    }
    if (test_rect_domain_contructor())
    {
        std::cout << "Test rectangular domain constructor succeeded!" << std::endl;
    }
    else
    {
        std::cout << "Test rectangular domain constuctor failed!" << std::endl;
    }
}






