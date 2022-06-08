#include "base_sweeper.h"
#include "../potential/harmonic_potential.h"
#include "../initial_condition/initial_condition.h"

class gSweeper:public BaseSweeper{
    public:
        gSweeper() = default;
        gSweeper(float start, float end, int num, bool endpoint, bool MPI_use, bool CUDA_use);
        void run(RectangularDomain *domain, InitialCondition *initial_condition, HarmonicPotential *potential);
    private: 
            
};