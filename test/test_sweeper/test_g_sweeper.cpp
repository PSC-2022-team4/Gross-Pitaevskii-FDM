#include "../../src/sweeper/g_sweeper.h"
#include "gtest/gtest.h"
#include <mpi.h>


TEST(GSweeperTest, gSweeperSerial){
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    if(rank == 0 ){
        bool all_passed = false;
        GSweeper gSweeper = GSweeper(float(0), float(-10), 2, false);
        ASSERT_FLOAT_EQ(gSweeper.get_value_from_idx(0), float(0));
        ASSERT_FLOAT_EQ(gSweeper.get_value_from_idx(1), float(-5));
        
        auto domain = RectangularDomain(11, 11, 0, 1e-8, 2, -1, 1, -1, 1);
        auto initial_cond_function = [](float x, float y)
        { return std::complex<float>{float(1.)}; }; //  expf(-(x * x + y * y) / (9))
        auto *initial_condition = new InitialCondition(initial_cond_function);
        auto *potential = new HarmonicPotential(3, 5);
        gSweeper.set_print_info(false);
        gSweeper.set_save_data(false);
        gSweeper.run(&domain,  initial_condition, potential);
        all_passed = true; 
        ASSERT_TRUE(all_passed);
        delete initial_condition;
        delete potential;
    }
    MPI_Barrier(comm);
}

TEST(GSweeperTest, gSweeperMPI){
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    bool passed = true;
    bool all_passed = false;
    
    GSweeper gSweeper = GSweeper(float(0), float(-10), 2, false);
    
    auto domain = RectangularDomain(21, 21, 0, 0.001, 2, -10, 10, -10, 10);
    auto initial_cond_function = [](float x, float y)
    { return std::complex<float>{float(1.)}; }; //  expf(-(x * x + y * y) / (9))
    
    auto *initial_condition = new InitialCondition(initial_cond_function);
    auto *potential = new HarmonicPotential(3, 5);
    gSweeper.set_MPI_info(rank, size);
    gSweeper.set_print_info(false);
    gSweeper.set_save_data(false);
    gSweeper.run(&domain,  initial_condition, potential);
    
    if (rank == 0)
    {
        int *passed_array = (int *)malloc(size * sizeof(int));
        //save results in passed_array
        MPI_Gather(&passed, 1, MPI_INT, passed_array, 1, MPI_INT, 0, MPI_COMM_WORLD);
        for (auto i = 0; i < size; ++i)
        {
            if (!passed_array[i])
            {

                break;
            }
            if (i == (size - 1))

            {
                all_passed = true;
            }
        }
    }
    else
    {
        MPI_Gather(&passed, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
    }
    ASSERT_TRUE(all_passed);
    delete initial_condition;
    delete potential;
}

TEST(GSweeperTest, gSweeperCUDASerial){
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    if(rank == 0 ){
        bool all_passed = false;
        GSweeper gSweeper = GSweeper(float(0), float(-10), 2, false);
        auto domain = RectangularDomain(21, 21, 0, 0.01, 2, -10, 10, -10, 10);
        auto initial_cond_function = [](float x, float y)
        { return std::complex<float>{float(1.)}; }; //  expf(-(x * x + y * y) / (9))
        auto *initial_condition = new InitialCondition(initial_cond_function);
        auto *potential = new HarmonicPotential(3, 5);
        gSweeper.set_print_info(false);
        gSweeper.set_save_data(false);
        gSweeper.set_CUDA_info(1,1);
        gSweeper.run(&domain,  initial_condition, potential);
        all_passed = true; 
        ASSERT_TRUE(all_passed);
        delete initial_condition;
        delete potential;
    }
    MPI_Barrier(comm);
    
}

TEST(GSweeperTest, gSweeperMPIandCUDA){
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    bool passed = true;
    bool all_passed = false;
    
    GSweeper gSweeper = GSweeper(float(0), float(-10), 2, false);
    
    auto domain = RectangularDomain(21, 21, 0, 0.001, 2, -10, 10, -10, 10);
    auto initial_cond_function = [](float x, float y)
    { return std::complex<float>{float(1.)}; }; //  expf(-(x * x + y * y) / (9))
    
    auto *initial_condition = new InitialCondition(initial_cond_function);
    auto *potential = new HarmonicPotential(3, 5);
    gSweeper.set_MPI_info(rank, size);
    gSweeper.set_CUDA_info(2, 3);
    gSweeper.set_print_info(false);
    gSweeper.set_save_data(false);
    gSweeper.run(&domain,  initial_condition, potential);
    
    if (rank == 0)
    {
        int *passed_array = (int *)malloc(size * sizeof(int));
        //save results in passed_array
        MPI_Gather(&passed, 1, MPI_INT, passed_array, 1, MPI_INT, 0, MPI_COMM_WORLD);
        for (auto i = 0; i < size; ++i)
        {
            if (!passed_array[i])
            {

                break;
            }
            if (i == (size - 1))

            {
                all_passed = true;
            }
        }
    }
    else
    {
        MPI_Gather(&passed, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
    }
    ASSERT_TRUE(all_passed);
    delete initial_condition;
    delete potential;
}
