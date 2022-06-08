#include "./base_sweeper.h"
#include "../potential/harmonic_potential.h"
#include "../initial_condition/initial_condition.h"

class GSweeper:public BaseSweeper{
    public:
        GSweeper() = default;
        GSweeper(float start, float end, int num, bool endpoint);
        void run(RectangularDomain *domain, InitialCondition *initial_condition, HarmonicPotential *potential);
    private: 
            
};