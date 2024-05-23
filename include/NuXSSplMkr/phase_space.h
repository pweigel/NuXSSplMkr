#ifndef __PHASE_SPACE_H
#define __PHASE_SPACE_H

#include "configuration.h"
#include "physconst.h"
#include "tools.h"
#include "math.h"

namespace nuxssplmkr {

class PhaseSpace {
    private:
        Configuration &config;
        
    public:
        PhaseSpace(Configuration &_cfg);
        ~PhaseSpace() {};
        void Initialize();
        void Compute(double m_fs);
        void SetFinalStateMass(double m_fs);

        bool Validate(double E);
        bool Validate(double E, double y);
        bool Validate(double E, double x, double y);

        void Print();

        PhysConst* pc; // Constants

        int flag = 0;
        bool debug = false;

        double Q2_min;
        double Q2_max;

        double W2_min;
        double W2_max;

        double x_min, x_max, y_min, y_max;
};

}
#endif