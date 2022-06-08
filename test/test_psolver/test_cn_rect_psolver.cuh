#pragma once
#include "../../src/solver/parallel_solver/crank_nicolson/cn_rect_psolver.cuh"

#include "../../src/potential/harmonic_potential.h"
#include "../../src/utils.h"
#include <functional>
#include <iostream>
#include <complex>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

void test_normalize(bool *passed);
void test_error_calculation(bool *passed);