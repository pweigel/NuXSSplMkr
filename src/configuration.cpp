#include "configuration.h"
using namespace std;

namespace nuxssplmkr{

Configuration::Configuration(string config_path) {
    // TODO: check if file exists
    // boost::property_tree::read_json(config_path, config);  // store settings in config
    std::ifstream infile(config_path);
    infile >> j;   
}

void Configuration::Populate() {
    unique_name = j.at("unique_name");
    pdfset = j.at("pdfset");
    replica = j.at("replica");

    mass_scheme = j.at("mass_scheme");

    int pto = j.at("perturbative_order"); // TODO: use existing enum
    perturbative_order = static_cast<QCDOrder>(pto);
    DIS_process = j.at("DIS_process");  // TODO: using existing enum
    current = CurrentMap.at(DIS_process);

    Nx    = j.at("Nx");
    NQ2   = j.at("NQ2");
    xmin  = j.at("xmin");
    xmax  = j.at("xmax");
    Q2min = j.at("Q2min");
    Q2max = j.at("Q2max");

    integral_min_Q2 = j.at("integral_Q2min");

    MassZ   = j.at("MassZ");
    MassW   = j.at("MassW");
    Rho     = j.at("Rho");
    Sin2ThW = j.at("Sin2ThW");

    Vud = j.at("Vud"); Vus = j.at("Vus"); Vub = j.at("Vub");
    Vcd = j.at("Vcd"); Vcs = j.at("Vcs"); Vcb = j.at("Vcb");
    Vtd = j.at("Vtd"); Vts = j.at("Vts"); Vtb = j.at("Vtb");

    // Make the pdf with LHAPDF and get its properties
    pdf = LHAPDF::mkPDF(pdfset, replica);
    for (int i=1; i<7; i++){
        pdf_quark_masses[i] = pdf->quarkMass(i);
    }
    PDFxmin  = pdf->xMin();
    PDFQ2min = pdf->q2Min();
    PDFQ2max = pdf->q2Max();

    if (current == CC) {
        M_boson2 = MassW * MassW;
    } else if (current == NC) {
        M_boson2 = MassZ * MassZ;
    }

    disable_top    = j.value("disable_top", false);
    enable_small_x = j.value("enable_small_x", false);
    small_x_order  = j.value("small_x_order", "NLL"); 
    evolve_pdf     = j.value("evolve_pdf", false);

}

void Configuration::Set_Target(string target_string) {
    target = target_string;
    target_type = TargetTypeMap.at(target_string);
}

void Configuration::Set_Projectile(string projectile_string) {
    projectile = projectile_string;
    neutrino_type = NeutrinoTypeMap.at(projectile_string);
    switch (neutrino_type) {
        case NeutrinoType::neutrino:
            cp_factor = 1.0;
            break;
        case NeutrinoType::antineutrino: 
            cp_factor = -1.0;
            break;
        default: 
            throw std::runtime_error("Cannot get CP Factor for specified projectile");
    }
}

void Configuration::Set_SF_Type(string _sf_type_string) {
    sf_type_string = _sf_type_string;
    sf_type = SFTypeMap.at(_sf_type_string);
}

void Configuration::LoadPDFSet() {
    
}

}