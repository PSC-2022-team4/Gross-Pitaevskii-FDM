#pragma once
#include <map>
#include <string>

struct DomainParameters
{

};
struct InitialConditionParameters
{
    std::string init_cond_type;
    std::map<std::string, float> init_cond_parameters;
};

struct EquationParameters
{
    float g;
    std::string potential_type;
    std::map<std::string, float> potential_parameters;
};

struct SolverParameters
{
};

struct Parameters
{
};