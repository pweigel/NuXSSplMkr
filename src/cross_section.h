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
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/foreach.hpp>
// #include <boost/math/interpolators/bilinear_uniform.hpp>
#include "interpolator.h"
#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/GridPDF.h"
// using boost::math::interpolators::bilinear_uniform;

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

struct Grid {
    int NQ2;
    int Nx;
    double Q2min, Q2max, xmin, xmax;
    std::vector<double> data;
};

class CrossSection {
    private:

        double APFEL_NC_Factor;
        double Normalization;
        double CP_Factor;

        bool F1_loaded = false;
        bool F2_loaded = false;
        bool F3_loaded = false;

        BilinearInterpolator F1_interpolator;
        BilinearInterpolator F2_interpolator;
        BilinearInterpolator F3_interpolator;

        double M_iso; // Isoscalar mass
        double M_l; // lepton mass
        double ENU; // neutrino energy


        double _kernel_y; // Used for integration

        // Limits of integration
        double integral_min_Q2;
        double integral_max_Q2;
        double integral_min_x;
        double integral_max_x;



        Configuration config;

    public:
        CrossSection(Configuration& _config);
        ~CrossSection() {};

        PhysConst* pc; // Constants

        // ~ Calculations ~
        double ds_dxdy_LO(double x, double y); // TODO: Remove or fix

        double ds_dxdy(double E, double x, double y);
        double ds_dxdy(double x, double y);
        double ds_dxdy_partonic(double E, double x, double y);
        double ds_dxdy_partonic(double x, double y);

        double ds_dy(double E, double y);
        double ds_dy_partonic(double E, double y);

        double ds_dxdy_TMC(); // TODO
        double ds_dy_TMC();   // TODO

        double TotalXS(double E);

        // ~ Kernel Functions ~
        double _ds_dy(double k);
        double ds_dxdy_kernel(double* k);
        double _ds_dy_partonic(double k);
        double _ds_dxdy_partonic(double* k);

        Grid Load_Grid(string path);

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