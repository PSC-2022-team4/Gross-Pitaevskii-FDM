#pragma once
#include "src/parallel_solver/crank_nicolson/cn_rect_psolver.cuh"
#include "src/utils.h"
#include <functional>
#include <iostream>
#include <complex>

void test_normalize(bool *passed);
void test_error_calculation(bool *passed);