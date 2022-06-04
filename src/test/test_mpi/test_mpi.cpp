#include "src/utils.h"
#include <mpi.h>
#include <iostream>

bool test_mpi_linking(int rank, int size){
    bool all_passed = false;
    MPI_Status status;
    int baton;
    if (rank == 0){
        std::cout<<"Total processor number:"<<size<<std::endl;
    }
    
    if (rank==0){
        baton = 1; 
        // baton 1개, MPI_INT 타입을 1번 프로세스에 tag = 999 로 전달 
        MPI_Send(&baton, 1, MPI_INT, 1, 999, MPI_COMM_WORLD);
        //baton 1개, MPI_INT 타입을 size-1 번 프로세스에서 받으려고 대기 
        MPI_Recv(&baton, 1, MPI_INT, size - 1, 999, MPI_COMM_WORLD, &status);
        all_passed = true;
    }else{
        // baton 받고 다음 프로세스에 던지기 
        MPI_Recv(&baton, 1, MPI_INT, rank-1, 999, MPI_COMM_WORLD, &status);
        MPI_Send(&baton, 1, MPI_INT, (rank+1)%size, 999, MPI_COMM_WORLD);
        all_passed = true;
    }
    return all_passed;
   
}

bool test_mpi_swap_g(int rank, int size){
    bool all_passed = true;
    
    double g = (double) rank / (double) size; 
    std::function<double(double, double)> potential;
    RectangularDomain *domain = (new RectangularDomain(21, 21, 0, 10, 11, -5, 5, -5, 5));
    auto initial_cond_function = [](double x, double y)
    { return std::complex<double>{1*std::exp(-(x*x + y*y)/(9))}; };

    auto *initial_condition =new  InitialCondition(initial_cond_function);
    initial_condition-> assign_to_domain(domain);
    potential= [](double x, double y ){
        return (double) 0.5 * (x*x + y *y);  };
    
    FERectSolver solver = FERectSolver(potential, g, domain);
    solver.solve();
    return all_passed;
}

void test_all_mpi(int rank, int size){
    bool passed = test_mpi_linking(rank, size);
    if(rank ==0){
        int *passed_array = (int*)malloc(size * sizeof(int));
        //save results in passed_array
        MPI_Gather(&passed, 1, MPI_INT, passed_array, 1, MPI_INT, 0, MPI_COMM_WORLD);
        for(auto i=0; i<size; ++i){
            if(!passed_array[i]){
                std::cout<<"MPI Test failed!!"<<std::endl;
                break;
            }
        }std::cout<<"MPI Test succeed"<<std::endl;

    }else{
        MPI_Gather(&passed, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
    }
    bool passed2 = test_mpi_swap_g(rank, size);
    if(rank ==0){
        int *passed_array = (int*)malloc(size * sizeof(int));
        //save results in passed_array
        MPI_Gather(&passed, 1, MPI_INT, passed_array, 1, MPI_INT, 0, MPI_COMM_WORLD);
        for(auto i=0; i<size; ++i){
            if(!passed_array[i]){
                std::cout<<"MPI Test with serial solver failed!!"<<std::endl;
                break;
            }
        }std::cout<<"MPI Test with serial solver succeed"<<std::endl;

    }else{
        MPI_Gather(&passed, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
    }
    
}