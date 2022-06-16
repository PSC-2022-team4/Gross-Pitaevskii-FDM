/**
 * @file config_parser.h
 * @author Gyeonghun Kim, Minyoung Kim
 * @brief Configuration Parser Class header.
 * @version 0.1
 * @date 2022-06-09
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#pragma once

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>

#include "parameters.h"

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
