#include <iostream>
#include "test/test_domain/test_base_domain.cpp"
#include "test/test_domain/test_rectangular_domain.cpp"
#include "test/test_armadillo/test_armadillo_linking.cpp"


int main(int argc, char **argv)
{
    test_all_base_domain();
    test_all_rectangular_domain();
    test_all_armadillo_linking();
    
}