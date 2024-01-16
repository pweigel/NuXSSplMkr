#ifndef __CONFIGURATION_H
#define __CONFIGURATION_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/filesystem.hpp>
#include "json.hpp"
#include "APFEL/APFEL.h"
#include "LHAPDF/LHAPDF.h"
#include <unordered_map>
using namespace std;

namespace nuxssplmkr {


// Enums
enum QCDOrder {LO, NLO, NNLO};
enum Current {CC, NC};
enum NeutrinoType {neutrino, antineutrino};
enum TargetType {proton, neutron};
enum Flavor {electron, muon, tau};
enum PDFVar {central, minus, plus};
enum SFType {total, light, charm, bottom, top};

// Enum maps
static unordered_map<string, QCDOrder> const QCDOrderMap = { {"LO",QCDOrder::LO}, {"NLO",QCDOrder::NLO}, {"NNLO",QCDOrder::NNLO} };
static unordered_map<string, Current> const CurrentMap = { {"CC",Current::CC}, {"NC",Current::NC} };
static unordered_map<string, NeutrinoType> const NeutrinoTypeMap = { {"neutrino",NeutrinoType::neutrino}, {"antineutrino",NeutrinoType::antineutrino} };
static unordered_map<string, TargetType> const TargetTypeMap = { {"proton",TargetType::proton}, {"neutron",TargetType::neutron} };
static unordered_map<string, Flavor> const FlavorMap = { {"electron",Flavor::electron},{"muon",Flavor::muon},{"tau",Flavor::tau} };
static unordered_map<NeutrinoType, double> CPFactorMap { {NeutrinoType::neutrino,1.},{NeutrinoType::antineutrino,-1.} };
static unordered_map<string, SFType> SFTypeMap { {"total", SFType::total},{"light", SFType::light},{"charm",SFType::charm},{"bottom",SFType::bottom},{"top",SFType::top} };

struct General_settings {
    string unique_name;
    bool debug;
};

struct PDF_settings {
    string pdfset;
    int replica;
    LHAPDF::PDF* pdf;
    std::map<int, double> pdf_quark_masses;
    double PDFxmin, PDFQ2min, PDFQ2max;
};

struct SF_settings {
    string mass_scheme;

    int pto;
    QCDOrder perturbative_order;
    string DIS_process;
    Current current;

    bool disable_top;
    bool enable_small_x;
    string small_x_order;
    bool evolve_pdf;
    bool enable_TMC;
    bool enable_CKMT;
    bool enable_PCAC;
    bool use_AlbrightJarlskog;
    
    int Nx, NQ2;
    double xmin, xmax, Q2min, Q2max;
};

struct CKMT_settings {
    double Q0;
    double Delta0;
    double AlphaR;
    double a, b, c, d;
    double F2A, F2B, F2f;
    double xF3A, xF3B, xF3f; // Note: we will apply cp_factor on xF3B!!!
};

struct xs_integration_settings {
    double xmin, xmax, Q2min, Q2max;
};

struct fundamental_constants {
  double MassZ, MassW, Rho, Sin2ThW;
  double Vud, Vus, Vub;
  double Vcd, Vcs, Vcb;
  double Vtd, Vts, Vtb;
  double Mboson2;
};

class Configuration {
  private:
    nlohmann::json j;

  public:
    // Configuration();
    Configuration(string config_path);
    ~Configuration() { };
    void Populate();
    void LoadPDFSet();
    void Set_Projectile(string projectile_string);
    void Set_Target(string target_string);
    void Set_SF_Type(string sf_type_string);

    // settings structs
    General_settings general;
    PDF_settings pdf;
    SF_settings SF;
    CKMT_settings CKMT;
    xs_integration_settings xs_integration;
    fundamental_constants constants;
    //

    string sf_type_string;
    SFType sf_type;
    string target;
    TargetType target_type;
    string projectile;
    NeutrinoType neutrino_type;

    double cp_factor;
};

}

#endif