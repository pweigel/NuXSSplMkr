#include "configuration.h"
#include "physconst.h"
#include "structure_function.h"
#include "cross_section.h"
#include "phase_space.h"

#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>

using namespace nuxssplmkr;

int main(int argc, char* argv[]){
    // if (argc < 5) {
    //     std::cout << "Not enough inputs!" << std::endl;
    //     return 1;
    // } else if (argc > 5) {
    //     std::cout << "Too many inputs!" << std::endl;
    //     return 1;
    // }
    const std::string config_path = argv[1]; // Path to .json file containing configuration info
    // const std::string projectile = argv[2]; // neutrino or antineutrino
    // const std::string target = argv[3]; // proton or neutron
    // const std::string sf_type = argv[4]; // Which SFs to use total, light, charm, ..

    std::cout << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "Config Path: " << config_path << std::endl;
    // std::cout << "Projectile: " << projectile << std::endl;
    // std::cout << "Target: " << target << std::endl;
    // std::cout << "SF Type: " << sf_type << std::endl;
    std::cout << "=============================================" << std::endl << std::endl;

    // double logemin = 2;
    // double logemax = 9;
    // int NE = 200;
    // double dE = (logemax - logemin) / (NE-1);

    PhysConst* pc = new PhysConst();

    Configuration config = Configuration(config_path);
    config.Populate();
    config.Set_Replica(config.pdf.replica);

    string sf_type = "charm";
    string projectile = "neutrino";
    string target = "proton";

    config.Set_SF_Type(sf_type);
    config.Set_Projectile(projectile);
    config.Set_Target(target);

    PhaseSpace ps(config);
    ps.Print();

    // std::cout << ps.Validate(100.0 * pc->GeV, 0.1, 0.1) << std::endl;
    // std::cout << "flag = " << ps.flag << std::endl;
    // std::cout << ps.Validate(100.0 * pc->GeV, 0.1, 1e-8) << std::endl;
    // std::cout << "flag = " << ps.flag << std::endl;

    std::string data_folder = "../data/" + config.general.unique_name + "/replica_" + std::to_string(config.pdf.replica);

    string f1 = data_folder + "/F1_" + projectile + "_" + target + "_" + sf_type + ".fits";
    string f2 = data_folder + "/F2_" + projectile + "_" + target + "_" + sf_type + ".fits";
    string f3 = data_folder + "/F3_" + projectile + "_" + target + "_" + sf_type + ".fits";

    CrossSection* xs = new CrossSection(config, ps);
    if (config.SF.mass_scheme != "parton") {
        xs->Load_Structure_Functions(f1, f2, f3);
    }
    xs->Set_Lepton_Mass(pc->muon_mass);

    double E, _xs;

    std::cout << "Checking the two different phase space methods..." << std::endl;

    E = 100.0 * pc->GeV;
    xs->use_phase_space = false;
    _xs = std::log10(xs->TotalXS(E));
    std::cout << "OLD (100 GeV): " << std::endl;
    std::cout << "  E [GeV] = " << E / pc->GeV << std::endl;
    std::cout << "  log10 sigma = " << _xs << std::endl; 
    xs->use_phase_space = true;
    _xs = std::log10(xs->TotalXS(E));
    std::cout << "NEW (100 GeV): " << std::endl;
    std::cout << "  E [GeV] = " << E / pc->GeV << std::endl;
    std::cout << "  log10 sigma = " << _xs << std::endl << std::endl; 

    E = 100.0 * pc->TeV;
    xs->use_phase_space = false;
    _xs = std::log10(xs->TotalXS(E));
    std::cout << "OLD (100 TeV): " << std::endl;
    std::cout << "  E [GeV] = " << E / pc->GeV << std::endl;
    std::cout << "  log10 sigma = " << _xs << std::endl; 
    xs->use_phase_space = true;
    _xs = std::log10(xs->TotalXS(E));
    std::cout << "NEW (100 TeV): " << std::endl;
    std::cout << "  E [GeV] = " << E / pc->GeV << std::endl;
    std::cout << "  log10 sigma = " << _xs << std::endl << std::endl; 

    std::cout << "Checking flags..." << std::endl;
    bool output;
    ps.debug = true;
    output = ps.Validate(100.0 * pc->GeV, 0.1, 0.1);
    std::cout << "  Flag 0: " << ps.flag << std::endl;
    output = ps.Validate(100.0 * pc->GeV, 0.1, 1e-5);
    std::cout << "  Flag 1: " << ps.flag << std::endl;
    output = ps.Validate(100.0 * pc->GeV, 1.1, 0.1);
    std::cout << "  Flag 2: " << ps.flag << std::endl;
    output = ps.Validate(100.0 * pc->GeV, 0.1, 1.1);
    std::cout << "  Flag 3: " << ps.flag << std::endl;
    config.Set_SF_Type("top");
    ps.Initialize();
    output = ps.Validate(100.0 * pc->GeV, 0.1, 0.1);
    std::cout << "  Flag 4: " << ps.flag << std::endl;


    return 0;
}