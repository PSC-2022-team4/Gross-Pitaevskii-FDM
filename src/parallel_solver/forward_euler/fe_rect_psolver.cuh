#pragma once
#include <complex>
#include <vector>
#include <functional>
#include "src/domain/rect_domain.h"
#include "src/initial_condition/initial_condition.h"
#include "src/serial_solver/base_serial_solver.h"
#include <cuda.h>
#include <cuda_runtime.h>

#define nTx 16
#define nTy 16

class FERectPSolver : public BaseSolver
{
public:
    FERectPSolver() = default;
    FERectPSolver(std::function<double(double, double)> potential,
                  double g,
                  RectangularDomain *domain);
    void solve(std::string dir_name="");
    void solve_single_time(int k);

protected:
    RectangularDomain *domain;
    RectangularSpatialGrid potential_grid;
    void generate_potential_grid();
    double get_potential_value(int i, int j);
    std::complex<double> temporal_equation(int i, int j, int k);
};

__global__ void fe_rect_cusolver(double *psi_old_real,
                                 double *psi_old_imag,
                                 double *psi_new_real,
                                 double *psi_new_imag,
                                 double *potential,
                                 int n_x,
                                 int n_y,
                                 double g,
                                 double h_x,
                                 double h_y,
                                 double tau);