#ifndef __STRUCTURE_FUNCTION_H
#define __STRUCTURE_FUNCTION_H
#include "configuration.h"
#include "physconst.h"
#include "LHAPDF/LHAPDF.h"

namespace nuxssplmkr {

class StructureFunction {
    private:
        Configuration config;

    public:
        // StructureFunction();
        StructureFunction(Configuration& _config);
        ~StructureFunction() {};

        PhysConst* pc; // Constants

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