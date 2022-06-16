/**
 * @file config_parser.h
 * @author Gyeonghun Kim, Minyoung Kim
 * @brief Configuration Parser Class header.
 * @version 0.1
 * @date 2022-06-9
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#pragma once

#include "parameters.h"
#include <string>

/**
 * @brief Configuration parser
 * 
 */
class ConfigParser{
    public:
        ConfigParser() = default;
        virtual ~ConfigParser() = default;
        static Parameters parse(std::string config_name, std::string filename);
        static Parameters get_default();
};
