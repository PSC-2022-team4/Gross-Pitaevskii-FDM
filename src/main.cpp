#include <iostream>
#include "test/test_domain/test_base_domain.cpp"

int main(int argc, char **argv)
{
    if(test_grid_point()){
        std::cout << "Test grid point constuctor succeeded!" << std::endl;
    }
    // if(test_base_domain_contructor()){
    //     std::cout << "Test base domain constructor succeeded!" << std::endl;
    // }
}