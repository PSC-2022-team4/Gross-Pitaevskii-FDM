/**
 * @file parameters.h
 * @author Gyeonghun Kim, Minyoung Kim
 * @brief Header for the Parameter classes.
 * @version 0.1
 * @date 2022-06-09
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#pragma once
#include <map>
#include <string>


/**
 * @brief Main branch parameters containing information about overall calculation
 * 
 */
class MainParameters
{
    public:
        std::string calculation_type = std::string("");
        std::map<std::string, float> float_parameters;
        std::map<std::string, int> int_parameters;
};

/**
 * @brief Domain branch parameters containing information about domain 
 * 
 */
class DomainParameters
{
    public:
    std::string domain_type = std::string("");
    int n_x, n_y, n_time;
    std::map<std::string, float> spatial_parameters;
    float time_start, time_end;
};

/**
 * @brief Initial condition branch parameters containing information about initial condition 
 * 
 */
class InitialConditionParameters
{
public:
    std::string init_cond_type = std::string("");
    std::map<std::string, float> init_cond_parameters;
};

/**
 * @brief Equation branch parameters containing information about equation 
 * 
 */
class EquationParameters
{
public:
    float g;
    std::string potential_type = std::string("");
    std::map<std::string, float> potential_parameters;
};

/**
 * @brief Solver branch parameters containing information about solver 
 * 
 */
class SolverParameters
{
public:
    std::string method = std::string("");
    bool run_parallel;
    bool save_data;
    bool print_info;
    std::map<std::string, float> solver_parameters;
    std::map<std::string, int> int_parameters;
};

/**
 * @brief Parameters containing configuration name and all other branched parameters
 * 
 */
class Parameters
{
public:
    std::string config_name;
    MainParameters main_parameters;
    DomainParameters domain_parameters;
    InitialConditionParameters init_cond_parameters;
    EquationParameters equation_parameters;
    SolverParameters solver_parameters;
};