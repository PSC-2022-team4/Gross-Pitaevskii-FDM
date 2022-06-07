#include "utils.h"

bool is_close(float a, float b, float tolerance){
    return std::abs(a - b) < tolerance;
}