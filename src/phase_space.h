#ifndef __PHASE_SPACE_H
#define __PHASE_SPACE_H

#include "configuration.h"

class PhaseSpace {
    private:
        double Q2_min;
        double Q2_max;

        double W2_min;
        double W2_max;
        
    public:
        PhaseSpace(Configuration _config);
        ~PhaseSpace();
        Compute();
        SetFinalStateMass();

        Configuration config;
        bool IsGood = false;  
}


#endif