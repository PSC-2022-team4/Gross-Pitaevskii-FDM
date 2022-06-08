#pragma once

#include "parameters.h"
#include <string>

class ConfigParser{
    public:
        ConfigParser() = default;
        virtual ~ConfigParser() = default;
        static Parameters parse(std::string config_name, std::string filename);
};
