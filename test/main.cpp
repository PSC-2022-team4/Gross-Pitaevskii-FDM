#include <mpi.h>
#include "gtest/gtest.h"
int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    ::testing::InitGoogleTest(&argc, argv);

    ::testing::TestEventListeners &listeners =
        ::testing::UnitTest::GetInstance()->listeners();
    if (rank != 0)
    {
        delete listeners.Release(listeners.default_result_printer());
    }
    auto result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}
