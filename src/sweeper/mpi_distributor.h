#include <mpi.h>

class MPIDistbutor{
    public:
        MPIDistbutor();
        MPIDistbutor(int rank, int size);
        void set_sweep_condition(float start, float end, int num, bool end_point);
        void check_run_condition();

    protected: 
        int rank;
        int size; 
        std::string information; 
        

}