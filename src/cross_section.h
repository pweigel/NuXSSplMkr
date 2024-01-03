#ifndef __CROSS_SECTION_H
#define __CROSS_SECTION_H
#include "configuration.h"
#include "structure_function.h"
#include "physconst.h"
#include "tools.h"
#include "math.h"
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_integration.h>
#include "photospline/splinetable.h"
#include <boost/math/quadrature/trapezoidal.hpp>
#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/GridPDF.h"

namespace nuxssplmkr {

template<class T,double (T::*f)(double*)>
double KernelHelper(double* x, size_t dim, void* param){
  T* p = (T*) param;
  return (p->*f)(x);
}

template<class T,double (T::*f)(double)>
double KernelHelper(double x, void* param){
  T* p = (T*) param;
  return (p->*f)(x);
}

template<class T,double (T::*f)(double,double),int n,int m>
double HK(double x, void* param){
  T* p = (T*) param;
  return ((double)m)*((p->*f)(x,p->Q2))/pow(x,(double)n);
}

class CrossSection {
    private:
        SFInfo sf_info;

        double APFEL_NC_Factor;
        double Normalization;
        double CP_Factor;

        bool F1_loaded = false;
        bool F2_loaded = false;
        bool F3_loaded = false;

        double M_iso; // Isoscalar mass
        double ENU; // neutrino energy

    public:
        CrossSection(Configuration config);
        ~CrossSection() {};

        PhysConst* pc; // Constants

        // ~ Calculations ~
        double ds_dxdy_LO(double x, double y);
        double ds_dxdy(double* k); // For integrators
        double ds_dxdy(double x, double y, double E);
        double ds_dxdy(double x, double y);
        double ds_dy();
        double ds_dxdy_TMC();
        double ds_dy_TMC();

        double TotalXS(double E);

        void Load_Structure_Functions(string sf1_path, string sf2_path, string sf3_path);
        void Load_F1(string path);
        void Load_F2(string path);
        void Load_F3(string path);

        photospline::splinetable<> F1;
        photospline::splinetable<> F2;
        photospline::splinetable<> F3;

        // ~ Settings ~
        void Set_Lepton_Mass(double m);
        void Set_Neutrino_Energy(double E);

};

}

#endif