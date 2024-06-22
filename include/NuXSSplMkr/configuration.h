#ifndef __CONFIGURATION_H
#define __CONFIGURATION_H

#include "NuXSSplMkr/physconst.h"
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
#include "LHAPDF/GridPDF.h"
#include <unordered_map>
using namespace std;

namespace nuxssplmkr {

// Enums
enum QCDOrder {LO, NLO, NNLO};
enum Current {CC, NC};
enum NeutrinoType {neutrino, antineutrino};
enum TargetType {proton, neutron, isoscalar};
enum Flavor {electron, muon, tau};
enum PDFVar {central, minus, plus};
enum SFType {total, light, charm, bottom, top};

// Enum maps
static unordered_map<string, QCDOrder> const QCDOrderMap = { {"LO",QCDOrder::LO}, {"NLO",QCDOrder::NLO}, {"NNLO",QCDOrder::NNLO} };
static unordered_map<string, Current> const CurrentMap = { {"CC",Current::CC}, {"NC",Current::NC} };
static unordered_map<string, NeutrinoType> const NeutrinoTypeMap = { {"neutrino",NeutrinoType::neutrino}, {"antineutrino",NeutrinoType::antineutrino} };
static unordered_map<string, TargetType> const TargetTypeMap = { {"proton",TargetType::proton}, {"neutron",TargetType::neutron}, {"isoscalar",TargetType::isoscalar}};
static unordered_map<string, Flavor> const FlavorMap = { {"electron",Flavor::electron},{"muon",Flavor::muon},{"tau",Flavor::tau} };
static unordered_map<NeutrinoType, double> CPFactorMap { {NeutrinoType::neutrino,1.},{NeutrinoType::antineutrino,-1.} };
static unordered_map<string, SFType> SFTypeMap { {"total", SFType::total},{"light", SFType::light},{"charm",SFType::charm},{"bottom",SFType::bottom},{"top",SFType::top} };

// codes
static unordered_map<string, string> SF_INTERACTION_CODE { {"CC", "1"}, {"NC", "2"} };
static unordered_map<string, string> SF_PARTICLE_CODE { {"neutrino", "0"}, {"antineutrino", "1"} };
static unordered_map<string, string> SF_NUMBER_CODES { {"F2", "1"}, {"FL", "2"}, {"xF3", "3'"}, {"F1", "4"}, {"F3", "5"} };
static unordered_map<string, string> SF_FLAVOR_CODES { {"total", "0"}, {"light", "1"}, {"charm", "2"}, {"bottom", "3"}, {"top", "4"} };

struct General_settings {
    string unique_name;
    string data_path;
    bool debug;
};

struct PDF_settings {
    string pdfset;
    int replica;
    std::map<int, double> pdf_quark_masses;
    double PDFxmin, PDFQ2min, PDFQ2max;
};

struct SF_settings {
    string mass_scheme;

    int pto;
    QCDOrder perturbative_order;
    string DIS_process;
    Current current;
    int FFNS;

    bool dynamic_fixed_flavor;
    bool enable_FONLL_damping;
    double FONLL_damping_factor;

    bool disable_top;
    int nf;
    bool enable_small_x;
    string small_x_order;
    bool evolve_pdf;
    bool enable_TMC;
    bool enable_CKMT;
    bool enable_PCAC;
    bool use_AlbrightJarlskog;
    double TMC_Q2max;
    
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

struct PCAC_settings {
    double A;
    double B;
};

struct xs_settings {
    bool enable_mass_terms;
    bool enable_shallow_region;
    bool enable_radiative_corrections;
    int mode;
    double xmin, xmax, ymin, ymax, Q2min, Q2max, W2min, W2max; // integration limits
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
    PhysConst* pc; // Constants
  
  public:
    // Configuration();
    Configuration(string config_path);
    ~Configuration() { };
    void Populate();
    void LoadPDFSet();
    void Set_Replica(int replica);
    void Set_Projectile(string projectile_string);
    void Set_Target(string target_string);
    void Set_SF_Type(string sf_type_string);
    void Set_Mass_Scheme(string mass_scheme);
    void Set_Lepton_Mass(double m);
    void Set_Mode(int mode);
    void Set_Perturbative_Order(int pto);
    string Get_SF_Code(string sf);
    LHAPDF::GridPDF* Get_LHAPDF_SF(int mode);

    // settings structs
    General_settings general;
    PDF_settings pdf_info;
    SF_settings SF;
    CKMT_settings CKMT;
    PCAC_settings PCAC;
    xs_settings XS;
    fundamental_constants constants;
    //

    LHAPDF::PDF* pdf;

    string sf_type_string;
    SFType sf_type;
    string target;
    TargetType target_type;
    string projectile;
    NeutrinoType neutrino_type;

    int mode;

    double cp_factor;
    double target_mass;
    double lepton_mass;


    bool flag_set_flavor = false;
    bool flag_set_projectile = false;
    bool flag_set_target = false;

    bool dynamic_W2_min = false;
};

}

#endif