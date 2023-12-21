#ifndef __CROSS_SECTION_H
#define __CROSS_SECTION_H
#include "configuration.h"
#include "structure_function.h"
#include "physconst.h"
#include "tools.h"
#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/GridPDF.h"
#include "APFEL/APFEL.h"
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_integration.h>
#include "photospline/splinetable.h"

namespace nuxssplmkr {

class CrossSection {
    private:
        StructureFunction sf;

        double APFEL_NC_Factor;
        double Normalization;
        double CP_Factor;

    public:
        CrossSection(StructureFunction& _sf);
        ~CrossSection() {};

        PhysConst* pc; // Constants

        // ~ Calculations ~
        double ds_dxdy(double x, double y, double Q2);
        double ds_dy();
        double ds_dxdy_TMC();
        double ds_dy_TMC();

        // double CrossSection(double x, double Q2);
        // double ds_dxdy(double x, double y, double Q2);
        // // double ds_dy(double x, double y, double Q2);
        // double Evaluate(double Q2, double x, double y);
        // double SigR_Nu_LO(double x, double y, map<int, double> dis);
        // double KernelXS(double * k);
        // double TotalXS();

        // ~ Settings ~
        void Set_Lepton_Mass(double m);
        void Set_Neutrino_Energy(double E);

};

}

#endif