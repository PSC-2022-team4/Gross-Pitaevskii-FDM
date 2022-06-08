#include <iostream>
#include <mpi.h>

#include "domain/rect_domain.h"
#include "initial_condition/initial_condition.h"
#include "potential/harmonic_potential.h"
#include "solver/base_solver.h"
#include "solver/parallel_solver/forward_euler/fe_rect_psolver.cuh"
#include "solver/parallel_solver/crank_nicolson/cn_rect_psolver.cuh"
#include "solver/serial_solver/forward_euler/fe_rect_solver.h"
#include "solver/serial_solver/crank_nicolson/cn_rect_solver.h"
#include "configuration/config_parser.h"
#include "configuration/parameters.h"
#include "sweeper/g_sweeper.h"
#include "sweeper/harmonic_p_sweeper.h"

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    if(argc != 2){
        std::cerr << "Input config file required. " << std::endl;
    }
    std::string config_filename = std::string(argv[1]);

    std::size_t dir_pos = config_filename.find_last_of("/");
    std::string dir = config_filename.substr(0, dir_pos);
    std::string config_name_extension = config_filename.substr(dir_pos+1, config_filename.length());
    std::string config_name = config_name_extension.substr(0, config_name_extension.length() - 4);

    auto parameters = ConfigParser::parse(config_name, config_filename);


    if (parameters.domain_parameters.domain_type == "rectangular")
    {
        if (parameters.main_parameters.calculation_type == "single")
        {
            RectangularDomain *domain = new RectangularDomain(parameters.domain_parameters.n_x,
                                                              parameters.domain_parameters.n_y,
                                                              parameters.domain_parameters.time_start,
                                                              parameters.domain_parameters.time_end,
                                                              parameters.domain_parameters.n_time,
                                                              parameters.domain_parameters.spatial_parameters["x_start"],
                                                              parameters.domain_parameters.spatial_parameters["x_end"],
                                                              parameters.domain_parameters.spatial_parameters["y_start"],
                                                              parameters.domain_parameters.spatial_parameters["y_end"]);

            std::function<std::complex<float>(float, float)> initial_cond_function;
            if (parameters.init_cond_parameters.init_cond_type == "singlegaussian")
            {
                auto sigma_x = parameters.init_cond_parameters.init_cond_parameters["sigma_x"];
                auto sigma_y = parameters.init_cond_parameters.init_cond_parameters["sigma_y"];
                auto center_x = parameters.init_cond_parameters.init_cond_parameters["center_x"];
                auto center_y = parameters.init_cond_parameters.init_cond_parameters["center_y"];
                initial_cond_function = [center_x, center_y, sigma_x, sigma_y](float x, float y)
                { 
                    return std::complex<float>{
                        float(1.) * expf(-((x - center_x) * (x - center_x) / (sigma_x * sigma_x) + (y - center_y) * (y - center_y) / (sigma_y * sigma_y)))
                        };
                };
                auto *initial_condition = new InitialCondition(initial_cond_function);
                initial_condition->assign_to_domain(domain);
            }
            else
            {
                std::cerr << "Unexpected initial condition" << std::endl;
            }
            if (parameters.equation_parameters.potential_type == "harmonic")
            {
                auto omega_x = parameters.equation_parameters.potential_parameters["omega_x"];
                auto omega_y = parameters.equation_parameters.potential_parameters["omega_y"];
                auto *potential = new HarmonicPotential(omega_x, omega_y);
                potential->calcualte_potential_in_grid(domain);
            }
            else
            {
                std::cerr << "Unexpected initial condition" << std::endl;
            }

            float g = parameters.equation_parameters.g;


            bool save_data = parameters.solver_parameters.save_data;
            bool print_info = parameters.solver_parameters.print_info;

            if (parameters.solver_parameters.method == "cranknicolson")
            {
                if(parameters.solver_parameters.run_parallel){
                    float converge_crit = parameters.solver_parameters.solver_parameters["converge_crit"];
                    int max_iter = parameters.solver_parameters.solver_parameters["max_iter"];
                    int cuda_device = parameters.solver_parameters.int_parameters["cuda_device"];
                    CNRectPSolver solver = CNRectPSolver(g, domain, cuda_device);
                    solver.solve(converge_crit, max_iter, std::to_string(rank), print_info, save_data);
                }
                else{
                    float converge_crit = parameters.solver_parameters.solver_parameters["converge_crit"];
                    int max_iter = parameters.solver_parameters.solver_parameters["max_iter"];
                    CNRectSolver solver = CNRectSolver(g, domain);
                    solver.solve(converge_crit, max_iter, std::to_string(rank), print_info, save_data);
                }
            }
            else if (parameters.solver_parameters.method == "forwardeuler"){
                if (parameters.solver_parameters.run_parallel)
                {
                    int cuda_device = parameters.solver_parameters.int_parameters["cuda_device"];
                    FERectPSolver solver = FERectPSolver(g, domain, cuda_device);
                    solver.solve(std::to_string(rank), print_info, save_data);
                }
                else
                {
                    FERectSolver solver = FERectSolver(g, domain);
                    solver.solve(std::to_string(rank), print_info, save_data);
                }
            }
        }
        else if (parameters.main_parameters.calculation_type == "g_sweep")
        {
            GSweeper gSweeper = GSweeper(parameters.main_parameters.float_parameters["sweep_start"],
                                         parameters.main_parameters.float_parameters["sweep_end"],
                                         parameters.main_parameters.int_parameters["sweep_count"],
                                         bool(parameters.main_parameters.int_parameters["end_point"]));
            RectangularDomain *domain = new RectangularDomain(parameters.domain_parameters.n_x,
                                                              parameters.domain_parameters.n_y,
                                                              parameters.domain_parameters.time_start,
                                                              parameters.domain_parameters.time_end,
                                                              parameters.domain_parameters.n_time,
                                                              parameters.domain_parameters.spatial_parameters["x_start"],
                                                              parameters.domain_parameters.spatial_parameters["x_end"],
                                                              parameters.domain_parameters.spatial_parameters["y_start"],
                                                              parameters.domain_parameters.spatial_parameters["y_end"]);

            std::function<std::complex<float>(float, float)> initial_cond_function;
            if (parameters.init_cond_parameters.init_cond_type == "singlegaussian")
            {
                auto sigma_x = parameters.init_cond_parameters.init_cond_parameters["sigma_x"];
                auto sigma_y = parameters.init_cond_parameters.init_cond_parameters["sigma_y"];
                auto center_x = parameters.init_cond_parameters.init_cond_parameters["center_x"];
                auto center_y = parameters.init_cond_parameters.init_cond_parameters["center_y"];
                initial_cond_function = [center_x, center_y, sigma_x, sigma_y](float x, float y)
                {
                    return std::complex<float>{
                        float(1.) * expf(-((x - center_x) * (x - center_x) / (sigma_x * sigma_x) + (y - center_y) * (y - center_y) / (sigma_y * sigma_y)))};
                };
                auto *initial_condition = new InitialCondition(initial_cond_function);
                initial_condition->assign_to_domain(domain);

                if (parameters.equation_parameters.potential_type == "harmonic")
                {
                    auto omega_x = parameters.equation_parameters.potential_parameters["omega_x"];
                    auto omega_y = parameters.equation_parameters.potential_parameters["omega_y"];
                    auto *potential = new HarmonicPotential(omega_x, omega_y);
                    potential->calcualte_potential_in_grid(domain);

                    bool save_data = parameters.solver_parameters.save_data;
                    bool print_info = parameters.solver_parameters.print_info;
                    bool cuda_use = (bool)parameters.main_parameters.int_parameters["cuda_use"];
                    bool mpi_use = (bool)parameters.main_parameters.int_parameters["mpi_use"];
                    if (cuda_use)
                    {
                        gSweeper.set_CUDA_info(parameters.main_parameters.int_parameters["gpu_count"],
                                            parameters.main_parameters.int_parameters["calculation_per_gpu"]);
                    }
                    if (mpi_use)
                    {
                        gSweeper.set_MPI_info(rank, size);
                    }
                    gSweeper.set_print_info(print_info);
                    gSweeper.set_save_data(save_data);

                    gSweeper.run(domain, initial_condition, potential);
                }
                else
                {
                    std::cerr << "Unexpected initial condition" << std::endl;
                }
            }
            else
            {
                std::cerr << "Unexpected initial condition" << std::endl;
            }
        }
        else if (parameters.main_parameters.calculation_type == "anisotropy_sweep")
        {

        }
    }
    else
    {
        std::cerr << "Unexpected Domain" << std::endl;
    }
    MPI_Finalize();
}
