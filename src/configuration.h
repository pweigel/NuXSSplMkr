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
#include <unordered_map>
using namespace std;

namespace nuxssplmkr {

// TODO: Find a better home for these
enum QCDOrder {LO, NLO, NNLO};
enum Current {CC, NC};
enum NeutrinoType {neutrino, antineutrino};
enum Flavor {electron, muon, tau};
enum PDFVar {central, minus, plus};

static unordered_map<string, QCDOrder> const QCDOrderMap = { {"LO",QCDOrder::LO}, {"NLO",QCDOrder::NLO}, {"NNLO",QCDOrder::NNLO} };
static unordered_map<string, Current> const CurrentMap = { {"CC",Current::CC}, {"NC",Current::NC} };
static unordered_map<string, NeutrinoType> const NeutrinoTypeMap = { {"neutrino",NeutrinoType::neutrino}, {"antineutrino",NeutrinoType::antineutrino} };
static unordered_map<string, Flavor> const FlavorMap = { {"electron",Flavor::electron},{"muon",Flavor::muon},{"tau",Flavor::tau} };
static unordered_map<NeutrinoType, double> CPFactorMap { {NeutrinoType::neutrino,1.},{NeutrinoType::antineutrino,-1} };

struct SFInfo {
    string pdfset;
    int replica;
    
    string mass_scheme;
    QCDOrder perturbative_order;
    string DIS_process;
    Current current;
    string projectile;
    NeutrinoType neutrino_type;
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