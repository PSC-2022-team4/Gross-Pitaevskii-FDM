#include "../../src/domain/base_domain.h"
#include "../../src/utils.h"
#include <iostream>
#include <experimental/filesystem>
#include <mpi.h>
#include "gtest/gtest.h"
namespace fs = std::experimental::filesystem;

TEST(BaseDomainTest, GridPointTest)
{
    auto grid_point = GridPoint(0.1, 0.1, std::complex<float>{10., 1.});
    ASSERT_FLOAT_EQ(grid_point.x, 0.1);
    ASSERT_FLOAT_EQ(grid_point.y, 0.1);
    ASSERT_FLOAT_EQ(grid_point.value.real(), 10.);
    ASSERT_FLOAT_EQ(grid_point.value.imag(), 1.);
}

TEST(BaseDomainTest, BaseSpatialGridTest)
{
    auto spatial_grid = BaseSpatialGrid(100, 100);
    ASSERT_FLOAT_EQ(spatial_grid.at(10, 10)->x, 0.);
    ASSERT_FLOAT_EQ(spatial_grid.at(10, 10)->y, 0.);
    ASSERT_FLOAT_EQ(spatial_grid.at(10, 10)->value.real(), 0.);
    ASSERT_FLOAT_EQ(spatial_grid.at(10, 10)->value.imag(), 0.);
}

TEST(BaseDomainTest, ConstructorTest)
{
    auto domain = BaseDomain(100, 100, 0., 10., 11);
    ASSERT_FLOAT_EQ(domain.get_t_start(), 0.);
    ASSERT_FLOAT_EQ(domain.get_t_end(), 10.);
    ASSERT_EQ(domain.get_num_times(), 11);
    ASSERT_FLOAT_EQ(domain.get_dt(), 1.);
}

TEST(BaseDomainTest, ExportFileTest)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    if (rank == 0){
        auto domain = BaseDomain(100, 100, 0., 10., 3);
        domain.generate_directory_name("test_initialize", false);
        domain.generate_single_txt_file(std::string("test_initialize"));
        int count = 0;
        for (const auto &file : fs::directory_iterator(domain.get_path()))
            count += 1;

        ASSERT_EQ(count, 1);
    }
}
