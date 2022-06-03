#include <iostream>
#include "test/test_domain/test_base_domain.cpp"
#include "test/test_domain/test_polar_domain.cpp"
#include "test/test_domain/test_rectangular_domain.cpp"
#include "test/test_initial_condition/test_initial_condition.cpp"

int main(int argc, char **argv)
{
    test_all_base_domain();
    test_all_rectangular_domain();
    test_all_polar_domain();
    test_initial_condition_rectangular();
    test_all_initial_condition();
}