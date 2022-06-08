#include "../../src/sweeper/base_sweeper.h"

#include "gtest/gtest.h"
#include <mpi.h>

TEST(BaseSweeperTest, initTest)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    if(rank == 0 ){
            
        BaseSweeper baseSweeper1 = BaseSweeper(float(0), float(10), 10 );
        BaseSweeper baseSweeper2 = BaseSweeper(float(0), float(9), 10, true);
        ASSERT_FLOAT_EQ(baseSweeper1.get_start(), float(0));
        ASSERT_FLOAT_EQ(baseSweeper1.get_end(), float(10));
        ASSERT_FLOAT_EQ(baseSweeper2.get_start(), float(0));
        ASSERT_FLOAT_EQ(baseSweeper2.get_end(), float(9));

        ASSERT_FLOAT_EQ(baseSweeper1.get_value_from_idx(9), float(9));
        ASSERT_FLOAT_EQ(baseSweeper2.get_value_from_idx(9), float(9));

        ASSERT_EQ(baseSweeper1.get_number_of_pts(), 10);
    }
    
};
