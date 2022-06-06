#include "utils.h"
#include <string>
#include <iostream>
#include <fstream>
bool is_close(double a, double b, double tolerance){
    return std::abs(a - b) < tolerance;
}

