#ifndef __STRUCTURE_FUNCTION_H
#define __STRUCTURE_FUNCTION_H

#include "configuration.h"
#include "physconst.h"
#include "tools.h"
#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/GridPDF.h"
#include "APFEL/APFEL.h"
#include <deque>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_integration.h>
#include <photospline/splinetable.h>

namespace nuxssplmkr {

template<class T,double (T::*f)(double,double),int n,int m>
double HK(double x, void* param){
  T* p = (T*) param;
  return ((double)m)*((p->*f)(x,p->_kernel_Q2))/pow(x,(double)n);
}

template<class T,double (T::*f)(double)>
double KernelWrapper(double x, void* param) {
    T* p = (T*) param;
    return (p->*f)(x);
}

class StructureFunction {
    private:
        double s_w, Lu2, Ld2, Ru2, Rd2;
        double M_lepton, d_lepton, d_nucleon;
        double M_iso, GF2;
        double ENU; // neutrino energy

        double CP_factor = std::numeric_limits<double>::max();

        Configuration config;

        bool has_LO_coefficients = false;
        bool splines_loaded = false;
        std::map<int,double> F2coef;
        std::map<int,double> F3coef;

        // PCAC constants (cite reno paper)
        double pion_mass = 135.0; // TODO: exact value
        double fpi = 0.93 * pion_mass;

        bool use_AlbrightJarlskog;

        double _Q2_cached = -1.0;

        photospline::splinetable<> spline_F1;
        photospline::splinetable<> spline_F2;
        photospline::splinetable<> spline_F3;

        int mode = 0;
        string insuffix = "";
        string outsuffix = "";

    public:
        StructureFunction(Configuration& _config);
        ~StructureFunction() {};

        PhysConst* pc; // Constants
        vector<int> partons {-6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 21};

        double _kernel_Q2;
        double _kernel_x;
        double _kernel_xi;

        string f1_grid_fn;
        string f2_grid_fn;
        string f3_grid_fn;

        // ~ APFEL stuff ~
        void InitializeAPFEL();

        // ~ Structure Function Calculatons ~
        void GetCoefficients();  // Get LO coefficients
        double F1(double x, double Q2);
        double F2(double x, double Q2);
        double FL(double x, double Q2);
        double F3(double x, double Q2);
        double xF3(double x, double Q2);
        double F4(double x, double Q2);
        double F5(double x, double Q2);

        double F2_LO(map<int, double>& xq_arr); // Calculate F2 from pdf
        double xF3_LO(map<int, double>& xq_arr); // Calculate xF3 from pdf

        double RescalingVariable(double Q2); // slow rescaling variable
        double NachtmannR(double x, double Q2);
        double NachtmannXi(double x, double Q2);
        double NachtmannXibar(double x, double Q2); // w/ mass rescaling

        template<class T,double (T::*f)(double,double),int n,int m>
        double HGeneric(double x, double Q2);

        double H2_kernel(double u);
        double H3_kernel(double u);
        double G2_kernel(double v);

        double H2(double xi, double Q2);
        double H3(double xi, double Q2);
        double G2(double xi, double Q2);

        double F1_TMC(double x, double Q2);
        double F2_TMC(double x, double Q2);
        double F3_TMC(double x, double Q2);
        double xF3_TMC(double x, double Q2);

        double CKMT_n(double Q2);
        double CKMT_Delta(double Q2);
        double F2_CKMT(double x, double Q2);
        double xF3_CKMT(double x, double Q2);
        double F3_CKMT(double x, double Q2);
        double F_CKMT(double x, double Q2, double A, double B, double f);
        double R_CKMT(double x, double Q2);
        double F1_CKMT(double _F2, double x, double Q2);

        double F1_PCAC(double x, double Q2);
        double FL_PCAC(double x, double Q2);
        double F2_PCAC(double x, double Q2);
        double F3_PCAC(double x, double Q2);

        void ConstructFONLL();
        std::tuple<double,double,double> EvaluateSFs(double x, double Q2);
        void Compute();

        void LoadGrids(string inpath);
        void LoadSplines(string inpath);
        
        void BuildGrids(string outpath);

        std::map<int,double> PDFExtract(double x, double Q2);
        
        // ~ Settings ~
        void Set_Mode(int mode);
        void Set_Target(string target_string);
        void Set_Projectile(string projectile_string);
        void Set_Lepton_Mass(double m);
        void Set_Neutrino_Energy(double E);
        void Set_SF_Type(string sf_type_string);
        void Set_Q_APFEL(double Q);
};

}

#endif