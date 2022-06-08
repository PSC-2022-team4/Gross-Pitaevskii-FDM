#include "../../src/domain/rect_domain.h"
#include "../../src/utils.h"
#include <iostream>
#include "gtest/gtest.h"

TEST(RectDomainTest, SpatialGridTest)
{
    auto spatial_grid = RectangularSpatialGrid(101, 101, 0, 10, 0, 20);

    ASSERT_FLOAT_EQ(spatial_grid.at(10, 10)->x, 1.);
    ASSERT_FLOAT_EQ(spatial_grid.at(10, 10)->y, 2.);
    ASSERT_FLOAT_EQ(spatial_grid.at(10, 10)->value.real(), 0.);
    ASSERT_FLOAT_EQ(spatial_grid.at(10, 10)->value.imag(), 0.);
}

TEST(RectDomainTest, RectDomainConstructorTest)
{
    auto domain = RectangularDomain(100, 100, 0., 10., 11, 0, 10, 0, 20);

    ASSERT_FLOAT_EQ(domain.get_t_start(), 0.);
    ASSERT_FLOAT_EQ(domain.get_t_end(), 10.);
    ASSERT_EQ(domain.get_num_times(), 11.);
    ASSERT_FLOAT_EQ(domain.get_dt(), 1.);
}
