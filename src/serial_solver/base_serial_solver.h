#pragma once
#include <complex>
#include <vector>
#include "src/domain/base_domain.h"

class InitialCondition{
    
};
class BaseSolver{
    public: 
        BaseSolver();
        BaseSolver(InitialCondition initialCondition);
    protected: 
        InitialCondition initialCondition;
        BaseDomain baseDomain; 



        
};