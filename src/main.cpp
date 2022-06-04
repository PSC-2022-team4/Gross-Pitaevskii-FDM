#include <iostream>

#include "src/test/test_domain/test_base_domain.cpp"
#include "src/test/test_domain/test_polar_domain.cpp"
#include "src/test/test_domain/test_rect_domain.cpp"
#include "src/test/test_initial_condition/test_initial_condition.cpp"
#include "src/test/test_armadillo/test_armadillo_linking.cpp"
#include "src/test/test_solver/test_base_serial_solver.cpp"
#include "test/test_solver/test_fe_rect_solver.cpp"
#include "src/test/test_solver/test_cn_rect_solver.cpp"
#include "src/test/test_psolver/test_cn_rect_psolver.cpp"
#include "src/test/test_psolver/test_fe_rect_psolver.cpp"
int main(int argc, char **argv)
{

    // test_all_base_domain();
    // test_all_rect_domain();
    // test_all_polar_domain();
    // test_all_initial_condition();
    // test_all_armadillo_linking();
    // test_all_base_serial_solver();
    // test_all_fe_rect_solver();
    // test_all_cn_rect_solver();
    // test_all_fe_rect_psolver();
    test_all_cn_rect_psolver();
}
