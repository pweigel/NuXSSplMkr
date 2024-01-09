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

// struct APFEL_settings {
//     string mass_scheme;
//     bool disable_top; // Used for CSMS calculation
//     bool pdf_evolution;
//     bool small_x_resummation;
//     string small_x_order;
// }

// struct sf_grid_settings {
//   int Nx, NQ2;
//   double xmin, xmax, Q2min, Q2max;
// }

// struct xs_integration_settings {
//     double xmin, xmax, Q2min, Q2max;
// }

// struct fundamental_constants {
//   double MassZ, MassW, Rho, Sin2ThW;
//   double Vud, Vus, Vub;
//   double Vcd, Vcs, Vcb;
//   double Vtd, Vts, Vtb;
// }

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

    string unique_name;
    string pdfset;
    int replica;
    
    string mass_scheme;
    QCDOrder perturbative_order;
    string DIS_process;
    Current current;
    string projectile;
    NeutrinoType neutrino_type;
    string target;
    TargetType target_type;

    SFType sf_type;
    string sf_type_string;

    int Nx, NQ2;
    double xmin, xmax, Q2min, Q2max;
    double MassZ, MassW, Rho, Sin2ThW;
    double Vud, Vus, Vub;
    double Vcd, Vcs, Vcb;
    double Vtd, Vts, Vtb;

    double integral_min_Q2;

    LHAPDF::PDF* pdf;
    std::map<int, double> pdf_quark_masses; 
    double PDFxmin;
    double PDFQ2min;
    double PDFQ2max;

    double M_boson2;  // squared mass of the boson (W or Z)
    double cp_factor;
    bool disable_top; // Set the top mass to ~bottom mass+0.1, used for CSMS calculation 
    bool enable_small_x;
    string small_x_order;
    bool evolve_pdf; // Evolve pdf from Q0 (typically mc, or 1.3 GeV)
    bool enable_pcac; // Use PCAC for low Q2
};

}

#endif