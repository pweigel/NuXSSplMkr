#ifndef __STRUCTURE_FUNCTION_H
#define __STRUCTURE_FUNCTION_H
#include "configuration.h"
#include "physconst.h"
#include "tools.h"
#include "LHAPDF/LHAPDF.h"
#include "APFEL/APFEL.h"

namespace nuxssplmkr {

class StructureFunction {
    private:
        double s_w, Lu2, Ld2, Ru2, Rd2;
        double Mw2, Mz2, M_iso, GF2;
        double M_boson2;
        Configuration config;

    public:
        // StructureFunction();
        StructureFunction(Configuration& _config);
        ~StructureFunction() {};

        PhysConst* pc; // Constants
        SFInfo sf_info;
        std::vector<int> parton_ids {-5,-4,-3,-2,-1,1,2,3,4,5,21};

        // ~ Apfel stuff ~
        void InitializeAPFEL();

        // ~ Structure Function Calculatons ~
        double F1(double x, double Q2);
        double F2(double x, double Q2);
        double xF3(double x, double Q2);
        double F3(double x, double Q2);
};

}

#endif