#include "../../src/sweeper/harmonic_p_sweeper.h"
#include "gtest/gtest.h"
#include <mpi.h>

TEST(HPSweeperTest, Serial){
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    if(rank == 0 ){
        bool all_passed = false;
        HPSweeper hpSweeper = HPSweeper(float(1), float(0.5), 2, true);
        ASSERT_FLOAT_EQ(hpSweeper.get_value_from_idx(0), float(1));
        ASSERT_FLOAT_EQ(hpSweeper.get_value_from_idx(1), float(0.5));
        
        auto domain = RectangularDomain(21, 21, 0, 0.01, 2, -10, 10, -10, 10);
        auto initial_cond_function = [](float x, float y)
        { return std::complex<float>{float(1.)}; }; //  expf(-(x * x + y * y) / (9))
        auto *initial_condition = new InitialCondition(initial_cond_function);
        float g = -1 ;
        hpSweeper.set_print_info(false);
        hpSweeper.set_save_data(false);
        hpSweeper.run(&domain, initial_condition, g);
        all_passed = true; 
        ASSERT_TRUE(all_passed);
        delete initial_condition;

    }
    MPI_Barrier(comm);
}

TEST(HPSweeperTest, MPI){
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    bool passed = true;
    bool all_passed = false;
    
    HPSweeper hpSweeper = HPSweeper(float(1), float(0.5), 2, true);
        
    auto domain = RectangularDomain(21, 21, 0, 0.001, 2, -10, 10, -10, 10);
    auto initial_cond_function = [](float x, float y)
    { return std::complex<float>{float(1.)}; }; //  expf(-(x * x + y * y) / (9))
    
    auto *initial_condition = new InitialCondition(initial_cond_function);
    float g = -1.;
    hpSweeper.set_MPI_info(rank, size);
    hpSweeper.set_print_info(false);
    hpSweeper.set_save_data(false);
    hpSweeper.run(&domain,  initial_condition, g);
        
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

}

TEST(HPSweeperTest, CUDASerial){
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    if(rank == 0 ){
        bool all_passed = false;
        HPSweeper hpSweeper = HPSweeper(float(1), float(0.5), 2, true);
        auto domain = RectangularDomain(21, 21, 0, 0.01, 2, -10, 10, -10, 10);
        auto initial_cond_function = [](float x, float y)
        { return std::complex<float>{float(1.)}; }; //  expf(-(x * x + y * y) / (9))
        auto *initial_condition = new InitialCondition(initial_cond_function);
        float g = -1; 
        hpSweeper.set_print_info(false);
        hpSweeper.set_save_data(false);
        hpSweeper.set_CUDA_info(1,1);
        hpSweeper.run(&domain,  initial_condition, g);
        all_passed = true; 
        ASSERT_TRUE(all_passed);
        delete initial_condition;

    }
    MPI_Barrier(comm);
    
}

TEST(HPSweeperTest, MPIandCUDA){
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    bool passed = true;
    bool all_passed = false;
    
    HPSweeper hpSweeper = HPSweeper(float(1), float(0.5), 2, true);
        
    auto domain = RectangularDomain(21, 21, 0, 0.001, 2, -10, 10, -10, 10);
    auto initial_cond_function = [](float x, float y)
    { return std::complex<float>{float(1.)}; }; //  expf(-(x * x + y * y) / (9))
    
    auto *initial_condition = new InitialCondition(initial_cond_function);
    float g = -1; 
    hpSweeper.set_CUDA_info(2, 3);
    hpSweeper.set_MPI_info(rank, size);
    hpSweeper.set_print_info(false);
    hpSweeper.set_save_data(false);
    hpSweeper.run(&domain,  initial_condition, g);
    
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
}
