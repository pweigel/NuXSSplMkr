#include "configuration.h"
using namespace std;

namespace nuxssplmkr{

Configuration::Configuration(string config_path) {
    // TODO: check if file exists
    boost::property_tree::read_json(config_path, config);  // store settings in config
}

void Configuration::Populate() {
    sf_info.pdfset = config.get<string>("pdfset");
    sf_info.replica = config.get<int>("replica");

    sf_info.mass_scheme = config.get<string>("mass_scheme");
    int pto = config.get<int>("perturbative_order"); // TODO: use existing enum
    sf_info.perturbative_order = static_cast<QCDOrder>(pto);
    sf_info.DIS_process = config.get<string>("DIS_process");  // TODO: using existing enum
    sf_info.current = CurrentMap.at(sf_info.DIS_process);
    
    // sf_info.projectile = config.get<string>("projectile");
    // sf_info.neutrino_type = NeutrinoTypeMap.at(sf_info.projectile);
    // sf_info.target = config.get<string>("target");
    // sf_info.target_type = TargetTypeMap.at(sf_info.target);

    sf_info.Nx = config.get<int>("Nx");
    sf_info.NQ2 = config.get<int>("NQ2");
    sf_info.xmin = config.get<double>("xmin");
    sf_info.xmax = config.get<double>("xmax");
    sf_info.Q2min = config.get<double>("Q2min");
    sf_info.Q2max = config.get<double>("Q2max");

    sf_info.MassZ = config.get<double>("MassZ");
    sf_info.MassW = config.get<double>("MassW");
    sf_info.Rho = config.get<double>("Rho");
    sf_info.Sin2ThW = config.get<double>("Sin2ThW");

    sf_info.Vud = config.get<double>("Vud"); sf_info.Vus = config.get<double>("Vus"); sf_info.Vub = config.get<double>("Vub");
    sf_info.Vcd = config.get<double>("Vcd"); sf_info.Vcs = config.get<double>("Vcs"); sf_info.Vcb = config.get<double>("Vcb");
    sf_info.Vtd = config.get<double>("Vtd"); sf_info.Vts = config.get<double>("Vts"); sf_info.Vtb = config.get<double>("Vtb");

    // Make the pdf with LHAPDF and get its properties
    sf_info.pdf = LHAPDF::mkPDF(sf_info.pdfset, sf_info.replica);
    for (int i=1; i<7; i++){
        sf_info.pdf_quark_masses[i] = sf_info.pdf->quarkMass(i);
    }
    sf_info.PDFxmin = sf_info.pdf->xMin();
    sf_info.PDFQ2min = sf_info.pdf->q2Min();
    sf_info.PDFQ2max = sf_info.pdf->q2Max();

    if (sf_info.current == CC) {
        sf_info.M_boson2 = sf_info.MassW * sf_info.MassW;
    } else if (sf_info.current == NC) {
        sf_info.M_boson2 = sf_info.MassZ * sf_info.MassZ;
    }

    // sf_info.sf_type_string = "total"; // default is total
    // sf_info.sf_type = SFType::total;
    sf_info.disable_top = config.get<bool>("disable_top");
    // sf_info.disable_top = false;
    sf_info.Use_APFEL_LO = true;  // TODO!
}

void Configuration::Set_Target(string target_string) {
    sf_info.target = target_string;
    sf_info.target_type = TargetTypeMap.at(target_string);
}

void Configuration::Set_Projectile(string projectile_string) {
    sf_info.projectile = projectile_string;
    sf_info.neutrino_type = NeutrinoTypeMap.at(projectile_string);
}

void Configuration::Set_SF_Type(string sf_type_string) {
    sf_info.sf_type_string = sf_type_string;
    sf_info.sf_type = SFTypeMap.at(sf_type_string);
}

void Configuration::LoadPDFSet() {
    
}

}