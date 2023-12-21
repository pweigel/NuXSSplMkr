#ifndef __STRUCTURE_FUNCTION_H
#define __STRUCTURE_FUNCTION_H
#include "configuration.h"
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

template<class T,double (T::*f)(double*)>
double KernelHelper(double* x,size_t dim, void* param){
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

class StructureFunction {
    private:
        double s_w, Lu2, Ld2, Ru2, Rd2;
        double Mw2, Mz2, M_iso, GF2;
        double M_boson2;
        double M_lepton, d_lepton, d_nucleon;
        double ENU; // neutrino energy

        double CP_factor = std::numeric_limits<double>::max();

        Configuration config;

        bool has_LO_coefficients = false;
        std::map<int,double> F2coef;
        std::map<int,double> F3coef;
        // Grids
        // std::vector<std::vector<double>> grid_F1;
        // std::vector<std::vector<double>> grid_F2;
        // std::vector<std::vector<double>> grid_xF3;

    public:
        StructureFunction(Configuration& _config);
        ~StructureFunction() {};

        PhysConst* pc; // Constants
        SFInfo sf_info;
        vector<int> partons {-5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21};

        // ~ APFEL stuff ~
        void InitializeAPFEL();

        // ~ Structure Function Calculatons ~
        void GetCoefficients();  // Get LO coefficients
        double F1(double x, double Q2);
        double F2(double x, double Q2);
        double FL(double x, double Q2);
        double F3(double x, double Q2);
        double xF3(double x, double Q2);

        double F2_LO(map<int, double>& xq_arr); // Calculate F2 from pdf
        // double F2_NLO(map<int, double>& xq_arr);
        double xF3_LO(map<int, double>& xq_arr); // Calculate xF3 from pdf

        void BuildGrids();

        std::map<int,double> PDFExtract(double x, double Q2);
        double CrossSection(double x, double Q2);
        double ds_dxdy(double x, double y, double Q2);
        // double ds_dy(double x, double y, double Q2);
        double Evaluate(double Q2, double x, double y);
        double SigR_Nu_LO(double x, double y, map<int, double> dis);
        double KernelXS(double * k);
        double TotalXS();

        // ~ Settings ~
        void Set_Lepton_Mass(double m);
        void Set_Neutrino_Energy(double E);
        void Set_Use_APFEL_LO(bool value);
        void Set_Q_APFEL(double Q);
};

}

#endif