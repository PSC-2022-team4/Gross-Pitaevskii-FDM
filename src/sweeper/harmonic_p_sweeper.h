/**
 * @file harmonic_p_sweeper.h
 * @author Minyoung Kim, Gyeonghun Kim
 * @brief Header for Harmonic potential strength sweeper class 
 * @version 0.1
 * @date 2022-06-10
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "./base_sweeper.h"
#include "../potential/harmonic_potential.h"
#include "../initial_condition/initial_condition.h"
/**
 * @brief In harmonic potential, sweep the strength of the potential in 1 dimension
 * 
 */
class HPSweeper:public BaseSweeper{
    public:
        HPSweeper() = default;
        HPSweeper(float start, float end, int num, bool endpoint);
        void run(RectangularDomain *domain, InitialCondition *initial_condition, float g );
    private: 
            
};