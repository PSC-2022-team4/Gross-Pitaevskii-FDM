#include "./base_sweeper.h"
#include "../potential/harmonic_potential.h"
#include "../initial_condition/initial_condition.h"

class HPSweeper:public BaseSweeper{
    public:
        HPSweeper() = default;
        HPSweeper(float start, float end, int num, bool endpoint);
        void run(RectangularDomain *domain, InitialCondition *initial_condition, float g );
    private: 
            
};