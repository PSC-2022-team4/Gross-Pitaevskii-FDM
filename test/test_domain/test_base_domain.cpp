#include "../../src/domain/base_domain.h"
#include "../../src/utils.h"
#include <iostream>
#include <experimental/filesystem>
#include "gtest/gtest.h"
namespace fs = std::experimental::filesystem;

TEST(DomainTest, GridPointTest)
{
    auto grid_point = GridPoint(0.1, 0.1, std::complex<float>{10., 1.});
    ASSERT_FLOAT_EQ(grid_point.x, 0.1);
    ASSERT_FLOAT_EQ(grid_point.y, 0.1);
    ASSERT_FLOAT_EQ(grid_point.value.real(), 10.);
    ASSERT_FLOAT_EQ(grid_point.value.imag(), 1.);
}

TEST(DomainTest, BaseSpatialGridTest)
{
    auto spatial_grid = BaseSpatialGrid(100, 100);
    ASSERT_FLOAT_EQ(spatial_grid.at(10, 10)->x, 0.);
    ASSERT_FLOAT_EQ(spatial_grid.at(10, 10)->y, 0.);
    ASSERT_FLOAT_EQ(spatial_grid.at(10, 10)->value.real(), 0.);
    ASSERT_FLOAT_EQ(spatial_grid.at(10, 10)->value.imag(), 0.);
}

TEST(DomainTest, BaseDomainConstructorTest)
{
    auto domain = BaseDomain(100, 100, 0., 10., 11);
    ASSERT_FLOAT_EQ(domain.get_t_start(), 0.);
    ASSERT_FLOAT_EQ(domain.get_t_end(), 10.);
    ASSERT_EQ(domain.get_num_times(), 11);
    ASSERT_FLOAT_EQ(domain.get_dt(), 1.);
}

TEST(DomainTest, BaseDomainExportFileTest)
{
    auto domain = BaseDomain(100, 100, 0., 10., 3);
    std::string directory_name = domain.generate_txt_file("test_initialize");
    int count = 0;
    for (const auto &file : fs::directory_iterator(directory_name))
        count += 1;

    ASSERT_EQ(count, 3);
}
