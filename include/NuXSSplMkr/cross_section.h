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

        int F1_code;
        int F2_code;
        int F3_code;

        LHAPDF::GridPDF* SF_PDF;

        BilinearInterpolator F1_interpolator;
        BilinearInterpolator F2_interpolator;
        BilinearInterpolator F3_interpolator;

        double M_iso; // Isoscalar mass
        double M_l; // lepton mass
        double ENU; // neutrino energy

        double kernel_y; // Used for integration
        double kernel_x; // Used for integration
        double kernel_L; // large log for rc integration

        // Limits of integration
        double integral_min_Q2;
        double integral_max_Q2;
        double integral_min_x;
        double integral_max_x;

        double min_W2;

        double charm_mass = 1.3;
        double bottom_mass = 4.75;
        double top_mass = 173.0; // for testing

        double rc_prefactor;

        Configuration &config;
        PhaseSpace &ps;

        bool rc_dsdxdy_loaded = false;

        void SetThresholdW2();

        int mode;
        string insuffix;
        string outsuffix;

    public:
        CrossSection(Configuration& _config, PhaseSpace& _ps);
        ~CrossSection() {};

        PhysConst* pc; // Constants

        // ~ Calculations ~
        double ds_dxdy_LO(double x, double y); // TODO: Remove or fix
        double ds_dxdQ2(double E, double x, double y);
        double ds_dxdQ2(double x, double y);
        double ds_dxdy(double E, double x, double y);
        double ds_dxdy(double x, double y);
        double ds_dxdy_partonic(double E, double x, double y);
        double ds_dxdy_partonic(double x, double y);

        double ds_dy(double E, double y);
        double ds_dy_partonic(double E, double y);

        double ds_dxdy_TMC(); // TODO
        double ds_dy_TMC();   // TODO

        double TotalXS(double E);
        double TotalXS_xQ2(double E);
        double AlternativeTotalXS(double E);
        std::tuple<double, double> dsdy_xlims(double s, double y);
        bool PhaseSpaceIsGood(double x, double y, double E);
        bool PhaseSpaceIsGood_Q2(double x, double Q2, double E);
        
        // ~ Kernel Functions ~
        double ds_dy_kernel(double k);
        double ds_dxdy_kernel(double* k);
        double ds_dxdQ2_kernel(double* k);
        double _ds_dy_partonic(double k);
        double _ds_dxdy_partonic(double* k);

        // ~ QED Corrections ~
        double qed_splitting(double z);  // splitting function
        double rc_jacobian(double x, double y, double z); // vars -> hat(vars)
        double rc_kernel(double k);
        double calculate_rc_dsdxdy(double z, double E, double x, double y);
        double calculate_rc_dsdxdy(double z, double E, double x, double y, double L);
        double rc_integrate(double E, double x, double y);
        void rc_load_dsdxdy(string spline_path); // load dsdxdy for rc
        photospline::splinetable<> rc_dsdxdy; // cross section for rc calcs TODO: rename this
        
        Grid Load_Grid(string path);

        void Load_Structure_Functions(string sf1_path, string sf2_path, string sf3_path);
        void Load_Structure_Functions(string inpath);
        void Load_F1(string path);
        void Load_F2(string path);
        void Load_F3(string path);

        photospline::splinetable<> F1;
        photospline::splinetable<> F2;
        photospline::splinetable<> F3;

        // FONLL complexity
        // 0: What APFEL gives
        // 1: "FONLL-C" for everything but t
        // 2: Full decomposition
        // int complexity = 0;
        // photospline::splinetable<> F1lZM, F1cZM, F1bZM, F1tZM;
        // photospline::splinetable<> F1cM, F1bM, F1tM;
        // photospline::splinetable<> F1cM0, F1bM0, F1tM0;

        // photospline::splinetable<> F2lZM, F2cZM, F2bZM, F2tZM;
        // photospline::splinetable<> F2cM, F2bM, F2tM;
        // photospline::splinetable<> F2cM0, F2bM0, F2tM0;

        // photospline::splinetable<> F3lZM, F3cZM, F3bZM, F3tZM;
        // photospline::splinetable<> F3cM, F3bM, F3tM;
        // photospline::splinetable<> F3cM0, F3bM0, F3tM0;

        // ~ Settings ~
        void Set_Mode(int _mode);
        void Set_Lepton_Mass(double m);
        void Set_Neutrino_Energy(double E);
};

}

#endif