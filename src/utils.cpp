#include "utils.h"
#include <string>
#include <iostream>
#include <fstream>
bool is_close(float a, float b, float tolerance){
    return std::abs(a - b) < tolerance;
}

