#include <iostream>
<<<<<<< HEAD
#include "src/test/test_domain/test_base_domain.cpp"
#include "src/test/test_domain/test_polar_domain.cpp"
#include "src/test/test_domain/test_rectangular_domain.cpp"
#include "src/test/test_initial_condition/test_initial_condition.cpp"
#include "src/test/test_armadillo/test_armadillo_linking.cpp"
// #include "test/test_solver/test_forward_euler_rectangular_solver.cpp"
#include "src/test/test_solver/test_base_serial_solver.cpp"
#include "src/test/test_solver/test_crank_nicolson_solver.cpp"
=======
#include "test/test_domain/test_base_domain.cpp"
#include "test/test_domain/test_polar_domain.cpp"
#include "test/test_domain/test_rectangular_domain.cpp"
#include "test/test_initial_condition/test_initial_condition.cpp"
#include "test/test_armadillo/test_armadillo_linking.cpp"
#include "test/test_solver/test_forward_euler_rectangular_solver.cpp"
#include "test/test_solver/test_base_serial_solver.cpp"
>>>>>>> ae9304929e6d3ef7fecba2d21225c249a7c3bd33

int main(int argc, char **argv)
{
    test_all_base_domain();
    test_all_rectangular_domain();
    test_all_polar_domain();
    test_initial_condition_rectangular();
    test_all_initial_condition();
    test_all_armadillo_linking();
    test_all_base_serial_solver();
<<<<<<< HEAD
    // test_forward_euler_rectangular_solver();
    // test_all_crank_nicolson_solver();
}

=======
    // test_all_forward_euler_rectangular_serial_solver();


}
>>>>>>> ae9304929e6d3ef7fecba2d21225c249a7c3bd33
