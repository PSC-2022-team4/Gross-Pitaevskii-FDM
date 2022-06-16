/**
 * @file g_sweeper.h
 * @author Minyoung Kim, Gyeonghun Kim
 * @brief Header for g parameter sweeper class 
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
 * @brief g parameter sweeper class that inherits Base sweeper class.
 * 
 */
class GSweeper:public BaseSweeper{
    public:
        GSweeper() = default;
        GSweeper(float start, float end, int num, bool endpoint);
        void run(RectangularDomain *domain, InitialCondition *initial_condition, HarmonicPotential *potential);
    private: 
            
};