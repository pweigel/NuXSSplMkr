#ifndef __CROSS_SECTION_H
#define __CROSS_SECTION_H
#include "configuration.h"
#include "structure_function.h"
#include "phase_space.h"
#include "physconst.h"
#include "tools.h"
#include "math.h"
#include <tuple>
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
#include <apfel/apfelxx.h>
#include "mlinterp.hpp"

namespace nuxssplmkr {

const std::vector<std::string> FLAVORS {
  "light", "charm", "bottom", "top"
};

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

// struct Grid {
//     int NQ2;
//     int Nx;
//     double Q2min, Q2max, xmin, xmax;
//     std::vector<double> data;
// };

class CrossSection {
    private:

        double APFEL_NC_Factor;
        double Normalization;
        double CP_Factor;

        // bool F1_loaded = false;
        // bool F2_loaded = false;
        // bool F3_loaded = false;

        int F1_code;
        int F2_code;
        int F3_code;

        LHAPDF::GridPDF* SF_PDF;

        double ENU; // neutrino energy

        double kernel_y; // Used for integration
        double kernel_x; // Used for integration
        double kernel_L; // large log for rc integration
        double kernel_born_dsdxdy; //used to speed up the rc calculation

        // Limits of integration
        double integral_min_Q2;
        double integral_max_Q2;
        double integral_min_x;
        double integral_max_x;

        double min_W2;

        double rc_prefactor;

        Configuration &config;
        PhaseSpace &ps;

        photospline::splinetable<> rc_spline;
        bool rc_spline_loaded = false;

        // Tools for interpolation
        int interp_nd[3] = {100, 120, 120};
        double interp_E[100];
        double interp_y[120];
        double interp_x[120];
        double interp_indata[1440000];
        bool interp_grid_loaded = false;

        void SetThresholdW2();

        int mode;
        string insuffix;
        string outsuffix;

    public:
        CrossSection(Configuration& _config, PhaseSpace& _ps);
        ~CrossSection() {};

        PhysConst* pc; // Constants

        // ~ Calculations ~
        double ds_dxdQ2(double E, double x, double y);
        double ds_dxdQ2(double x, double y);
        double ds_dxdy(double E, double x, double y);
        double ds_dxdy(double x, double y);

        double ds_dy(double E, double y);

        double TotalXS(double E);
        double TotalXS_xQ2(double E);
        bool PhaseSpaceIsGood(double x, double y, double E);
        bool PhaseSpaceIsGood_Q2(double x, double Q2, double E);
        
        // ~ Kernel Functions ~
        double ds_dy_kernel(double k);
        double ds_dxdy_kernel(double* k);
        double ds_dxdQ2_kernel(double* k);

        // ~ QED Corrections ~
        double qed_splitting(double z);  // splitting function
        double qed_splitting_integrated(double a); // splitting function integrated from 0 to a
        double rc_jacobian(double x, double y, double z); // vars -> hat(vars)
        double rc_kernel(double k);
        double calculate_rc_dsdzdxdy(double z, double x, double y);
        double rc_dsdxdy(double E, double x, double y, double born_dsdxdy);
        void Load_RC_Spline(string spline_path);
        void Load_InterpGrid(string grid_path);
    
        // ~ Settings ~
        void Set_Mode(int _mode);
        void Set_Neutrino_Energy(double E);
};

}

#endif