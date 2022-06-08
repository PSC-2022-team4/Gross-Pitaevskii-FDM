#include "../../src/utils.h"
#include "../../src/potential/harmonic_potential.h"
#include "../../src/solver/serial_solver/forward_euler/fe_rect_solver.h"

#include <mpi.h>
#include <iostream>
#include "gtest/gtest.h"

TEST(MPITest, Linking)
{
    bool all_passed = false;
    MPI_Status status;
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    int baton;
    if (rank == 0)
    {
        std::cout << "Total processor number:" << size << std::endl;
    }

    if (rank == 0)
    {
        baton = 1;
        // baton 1개, MPI_INT 타입을 1번 프로세스에 tag = 999 로 전달
        MPI_Send(&baton, 1, MPI_INT, 1, 999, MPI_COMM_WORLD);
        //baton 1개, MPI_INT 타입을 size-1 번 프로세스에서 받으려고 대기
        MPI_Recv(&baton, 1, MPI_INT, size - 1, 999, MPI_COMM_WORLD, &status);
        all_passed = true;
    }
    else
    {
        // baton 받고 다음 프로세스에 던지기
        MPI_Recv(&baton, 1, MPI_INT, rank - 1, 999, MPI_COMM_WORLD, &status);
        MPI_Send(&baton, 1, MPI_INT, (rank + 1) % size, 999, MPI_COMM_WORLD);
        // all_passed = true;
    }

    if (rank == 0)
    {
        int *passed_array = (int *)malloc(size * sizeof(int));
        //save results in passed_array
        MPI_Gather(&all_passed, 1, MPI_INT, passed_array, 1, MPI_INT, 0, MPI_COMM_WORLD);
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
        MPI_Gather(&all_passed, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
    }
    ASSERT_TRUE(all_passed);
}

TEST(MPITest, SwapG)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    bool passed = true;
    bool all_passed = false;

    float g = (float)rank;

    RectangularDomain *domain = (new RectangularDomain(21, 21, 0, 5, 101, -5, 5, -5, 5));
    auto initial_cond_function = [](float x, float y)
    { return std::complex<float>{1 * std::exp(-(x * x + y * y) / (9))}; };

    auto *initial_condition = new InitialCondition(initial_cond_function);
    initial_condition->assign_to_domain(domain);
    auto *potential = new HarmonicPotential(3, 5);
    potential->calcualte_potential_in_grid(domain);

    FERectSolver solver = FERectSolver(g, domain);
    solver.solve(std::to_string(rank), false, false);

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
}
