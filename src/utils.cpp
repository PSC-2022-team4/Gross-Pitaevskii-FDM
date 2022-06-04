#include "utils.h"

bool is_close(double a, double b, double tolerance){
    return std::abs(a - b) < tolerance;
}