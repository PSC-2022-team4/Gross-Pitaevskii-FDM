#include <armadillo>
#include "../../src/utils.h"

bool test_armadillo_matrix_initialize()
{
    //Make 2 by 3 zeros matrix
    arma::mat A(2, 3, arma::fill::zeros);
    bool all_passed = true;
    if (A.n_rows != 2)
    {
        all_passed = false;
    }
    if (A.n_cols != 3)
    {
        all_passed = false;
    }
    if (!is_close(A(0, 0), 0., 1e-12))
    {
        all_passed = false;
    }
    if (!is_close(A(0, 1), 0., 1e-12))
    {
        all_passed = false;
    }
    if (!is_close(A(1, 0), 0., 1e-12))
    {
        all_passed = false;
    }
    if (!is_close(A(1, 1), 0., 1e-12))
    {
        all_passed = false;
    }
    return all_passed;
}
// bool test_armadillo_matrix(){

// }
bool test_all_armadillo_linking()
{
    if (test_armadillo_matrix_initialize())
    {
        std::cout << "Test armadillo matrix initialization succeeded!" << std::endl;
    }
    else
    {
        std::cout << "Test armadillo matrix initialization failed!" << std::endl;
    }
}