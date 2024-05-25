#include "NuXSSplMkr/configuration.h"
using namespace std;

namespace nuxssplmkr{

Configuration::Configuration(string config_path) {
    // TODO: check if file exists
    // boost::property_tree::read_json(config_path, config);  // store settings in config
    std::ifstream infile(config_path);
    infile >> j;   
}

void Configuration::Populate() {
    general.unique_name = j["general"].at("unique_name");
    general.data_path = j["general"].value("data_path", "../data");
    general.debug = j["general"].value("debug", false);

    SF.mass_scheme = j["SF"].at("mass_scheme");
    SF.pto = j["SF"].at("perturbative_order");
    SF.perturbative_order = static_cast<QCDOrder>(SF.pto);
    SF.DIS_process = j["SF"].at("DIS_process");
    SF.current = CurrentMap.at(SF.DIS_process);
    SF.FFNS = -1; // TODO: Figure out how I want to handle this

    SF.enable_FONLL_damping = j["SF"].value("enable_FONLL_damping", true);
    SF.FONLL_damping_factor = j["SF"].value("FONLL_damping_factor", 2.0); // APFEL default is 2


    SF.disable_top = j["SF"].value("disable_top", false);
    SF.enable_small_x = j["SF"].value("enable_small_x", false);
    SF.small_x_order = j["SF"].value("small_x_order", "NLL");
    SF.evolve_pdf = j["SF"].value("evolve_pdf", false);
    SF.enable_TMC = j["SF"].value("enable_TMC", false);
    SF.enable_CKMT = j["SF"].value("enable_CKMT", false);
    SF.enable_PCAC = j["SF"].value("enable_PCAC", false);
    SF.TMC_Q2max = j["SF"].value("TMC_Q2max", 30.0);
    SF.use_AlbrightJarlskog = j["SF"].value("use_AlbrightJarlskog", true);

    SF.Nx = j["SF"].at("Nx");
    SF.NQ2 = j["SF"].at("NQ2");
    SF.xmin = j["SF"].at("xmin");
    SF.xmax = j["SF"].at("xmax");
    SF.Q2min = j["SF"].at("Q2min");
    SF.Q2max = j["SF"].at("Q2max");

    CKMT.Q0 = j["CKMT"].value("Q0", 2.0); // Matching scale
    CKMT.Delta0 = j["CKMT"].value("Delta0", 0.07684);
    CKMT.AlphaR = j["CKMT"].value("AlphaR", 0.4150);
    CKMT.a = j["CKMT"].value("a", 0.2631);
    CKMT.b = j["CKMT"].value("b", 0.6452);
    CKMT.c = j["CKMT"].value("c", 3.5489);
    CKMT.d = j["CKMT"].value("d", 1.1170);
    CKMT.F2A = j["CKMT"]["F2"].value("A", 0.5967);   CKMT.F2B  = j["CKMT"]["F2"].value("B", 2.7145);   CKMT.F2f  = j["CKMT"]["F2"].value("f", 0.5962);
    CKMT.xF3A = j["CKMT"]["xF3"].value("A", 9.3955e-3); CKMT.xF3B = j["CKMT"]["xF3"].value("B", 2.4677);  CKMT.xF3f = j["CKMT"]["xF3"].value("f", 0.5962);

    PCAC.A = j["PCAC"].value("A", 0.147);
    PCAC.B = j["PCAC"].value("B", 0.265);

    XS.mode = j["XS"].value("mode", 1);
    XS.enable_mass_terms = j["XS"].value("enable_mass_terms", false);
    XS.enable_shallow_region = j["XS"].value("enable_shallow_region", false);
    XS.xmin = j["XS"]["integration"].at("xmin");
    XS.xmax = j["XS"]["integration"].at("xmax");
    XS.Q2min = j["XS"]["integration"].at("Q2min");
    XS.Q2max = j["XS"]["integration"].at("Q2max");

    constants.MassZ = j["constants"].at("MassZ");
    constants.MassW = j["constants"].at("MassW");
    constants.Rho = j["constants"].at("Rho");
    constants.Sin2ThW = j["constants"].at("Sin2ThW");
    constants.Vud = j["constants"].at("Vud"); constants.Vus = j["constants"].at("Vus"); constants.Vub = j["constants"].at("Vub");
    constants.Vcd = j["constants"].at("Vcd"); constants.Vcs = j["constants"].at("Vcs"); constants.Vcb = j["constants"].at("Vcb");
    constants.Vtd = j["constants"].at("Vtd"); constants.Vts = j["constants"].at("Vts"); constants.Vtb = j["constants"].at("Vtb");

    // Make the pdf with LHAPDF and get its properties
    pdf.pdfset = j["PDF"].at("pdfset");
    pdf.replica = j["PDF"].at("replica");
    pdf.pdf = LHAPDF::mkPDF(pdf.pdfset, pdf.replica);
    for (int i=1; i<7; i++){
        pdf.pdf_quark_masses[i] = pdf.pdf->quarkMass(i);
    }
    pdf.PDFxmin  = pdf.pdf->xMin();
    pdf.PDFQ2min = pdf.pdf->q2Min();
    pdf.PDFQ2max = pdf.pdf->q2Max();

    if (SF.current == CC) {
        constants.Mboson2 = constants.MassW * constants.MassW;
    } else if (SF.current == NC) {
        constants.Mboson2 = constants.MassZ * constants.MassZ;
    }

    dynamic_W2_min = j["XS"].value("dynamic_W2_min", false);
}

void Configuration::Set_Replica(int replica) {
    // Make the pdf with LHAPDF and get its properties
    pdf.replica = replica;
    pdf.pdf = LHAPDF::mkPDF(pdf.pdfset, pdf.replica);
    for (int i=1; i<7; i++){
        pdf.pdf_quark_masses[i] = pdf.pdf->quarkMass(i);
    }
    pdf.PDFxmin  = pdf.pdf->xMin();
    pdf.PDFQ2min = pdf.pdf->q2Min();
    pdf.PDFQ2max = pdf.pdf->q2Max();
}

void Configuration::Set_Target(string target_string) {
    target = target_string;
    target_type = TargetTypeMap.at(target_string);

    flag_set_target = true;
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
    flag_set_projectile = true;
}

void Configuration::Set_SF_Type(string _sf_type_string) {
    sf_type_string = _sf_type_string;
    sf_type = SFTypeMap.at(_sf_type_string);
    // If we use FONLL, we want it to use the right threshold behavior
    // For example, for Fc we want to do M3 -> M3 - D(ZM4 - M03) at m_c^2
    if (SF.mass_scheme == "FONLL-A" || SF.mass_scheme == "FONLL-B" || SF.mass_scheme == "FONLL-C") {
        switch (sf_type) {
            case SFType::charm : SF.FFNS = 3 ; break;
            case SFType::bottom: SF.FFNS = 4 ; break;
            case SFType::top   : SF.FFNS = 5 ; break;
            default            : SF.FFNS = -1; break;
        }
        if (SF.FFNS > 0) {
            std::cout << "WARNING: FFNS is set to " << SF.FFNS << " to properly match FONLL calculation!" << std::endl;
        }
    }
    flag_set_flavor = true;
}

void Configuration::Set_Mass_Scheme(string mass_scheme) {
    SF.mass_scheme = mass_scheme;
}

void Configuration::Set_Perturbative_Order(int pto) {
    SF.pto = pto;
    SF.perturbative_order = static_cast<QCDOrder>(SF.pto);
}

void Configuration::LoadPDFSet() {
    
}

}