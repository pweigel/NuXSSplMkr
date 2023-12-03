#ifndef __CONFIGURATION_H
#define __CONFIGURATION_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/filesystem.hpp>
#include "APFEL/APFEL.h"
#include "LHAPDF/LHAPDF.h"
using namespace std;

namespace nuxssplmkr {

struct SFInfo {
    string pdfset;
    int replica;
    
    string mass_scheme;
    int perturbative_order;
    string DIS_process;
    string projectile;
    string target;

    int Nx, NQ2;
    double xmin, xmax, Q2min, Q2max;
    double MassZ, MassW, Rho, Sin2ThW;
    double Vud, Vus, Vub;
    double Vcd, Vcs, Vcb;
    double Vtd, Vts, Vtb;

    LHAPDF::PDF* pdf;
    std::map<int, double> pdf_quark_masses; 
    double PDFxmin;
    double PDFQ2min;
    double PDFQ2max;
};

class Configuration {
  private:
    boost::property_tree::ptree config;

  public:
    // Configuration();
    Configuration(string config_path);
    ~Configuration() { };
    void Populate();
    void LoadPDFSet();
    SFInfo sf_info;
};

}

#endif