/**
 * @file rect_domain.h
 * @author Minyoung Kim, Gyeonghun Kim
 * @brief Header for Rectangular domain class
 * @version 0.1
 * @date 2022-06-02
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#pragma once
#include "base_domain.h"

/**
 * @brief Rectangular Spatial Grid class for single time step.
 * 
 */
class RectangularSpatialGrid : public BaseSpatialGrid
{
public:
    RectangularSpatialGrid() = default;
    RectangularSpatialGrid(int num_grid_1, int num_grid_2, float x_start, float x_end, float y_start, float y_end);
    ~RectangularSpatialGrid();

private:
    float x_start;
    float x_end;
    float y_start;
    float y_end;
};

/**
 * @brief Rectangular domain containing multiple timesteps.
 * 
 */
class RectangularDomain : public BaseDomain
{
public:
    RectangularDomain() = default;
    RectangularDomain(int num_grid_1, int num_grid_2, float t_start, float t_end, int num_times, float x_start, float x_end, float y_start, float y_end);
    float get_x_start();
    float get_y_start();
    float get_x_end();
    float get_y_end();
    void generate_single_txt_file(std::string filename, bool cuda_mode = false);
    // a grid to save potential values
    RectangularSpatialGrid *potential_grid;
    void update_time(bool cuda_mode = false);
    void reset();
    ~RectangularDomain();

private:
    float x_start;
    float x_end;
    float y_start;
    float y_end;
};