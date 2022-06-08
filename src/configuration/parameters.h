#pragma once
#include <map>
#include <string>

class MainParameters
{
    public:
        std::string calculation_type = std::string("");
        std::map<std::string, float> float_parameters;
        std::map<std::string, int> int_parameters;
};
class DomainParameters
{
    public:
    std::string domain_type = std::string("");
    int n_x, n_y, n_time;
    std::map<std::string, float> spatial_parameters;
    float time_start, time_end;
};
class InitialConditionParameters
{
public:
    std::string init_cond_type = std::string("");
    std::map<std::string, float> init_cond_parameters;
};

class EquationParameters
{
public:
    float g;
    std::string potential_type = std::string("");
    std::map<std::string, float> potential_parameters;
};

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