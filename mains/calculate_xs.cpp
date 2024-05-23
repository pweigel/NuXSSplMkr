#include "NuXSSplMkr/configuration.h"
#include "NuXSSplMkr/physconst.h"
#include "NuXSSplMkr/phase_space.h"
#include "NuXSSplMkr/structure_function.h"
#include "NuXSSplMkr/cross_section.h"
#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>

using namespace nuxssplmkr;

int main(int argc, char* argv[]){
    if (argc < 5) {
        std::cout << "Not enough inputs!" << std::endl;
        return 1;
    } else if (argc > 5) {
        std::cout << "Too many inputs!" << std::endl;
        return 1;
    }
    const std::string config_path = argv[1]; // Path to .json file containing configuration info
    const std::string projectile = argv[2]; // neutrino or antineutrino
    const std::string target = argv[3]; // proton or neutron
    const std::string sf_type = argv[4]; // Which SFs to use total, light, charm, ..

    std::cout << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "Config Path: " << config_path << std::endl;
    std::cout << "Projectile: " << projectile << std::endl;
    std::cout << "Target: " << target << std::endl;
    std::cout << "SF Type: " << sf_type << std::endl;
    std::cout << "=============================================" << std::endl << std::endl;

    double logemin = 1;
    double logemax = 9;
    int NE = 200;
    double dE = (logemax - logemin) / (NE-1);

    PhysConst* pc = new PhysConst();

    Configuration config = Configuration(config_path);
    config.Populate();
    config.Set_Replica(config.pdf.replica);
    std::string data_folder = "../data/" + config.general.unique_name + "/replica_" + std::to_string(config.pdf.replica);
    std::cout << "Loading/saving data to: " << data_folder << std::endl;

    // Make the cross sections folder if it doesn't exist
    boost::filesystem::path out_folder = data_folder + "/cross_sections/";
    if (!boost::filesystem::exists(out_folder)) {
        boost::filesystem::create_directories(out_folder);
    }
    
    config.Set_SF_Type(sf_type);
    std::cout << "SFType set to: " << sf_type << std::endl;

    config.Set_Projectile(projectile);
    std::cout << "Projectile set to: " << projectile << std::endl;

    config.Set_Target(target);
    std::cout << "Target set to: " << target << std::endl;

    PhaseSpace ps(config);
    ps.Print();

    // string f1 = data_folder + "/F1_" + projectile + "_" + target + "_" + sf_type + ".fits";
    // string f2 = data_folder + "/F2_" + projectile + "_" + target + "_" + sf_type + ".fits";
    // string f3 = data_folder + "/F3_" + projectile + "_" + target + "_" + sf_type + ".fits";

    CrossSection* xs = new CrossSection(config, ps);
    if (config.SF.mass_scheme != "parton") {
        xs->Load_Structure_Functions(data_folder);
    }
    xs->Set_Lepton_Mass(pc->muon_mass);

    // total
    std::ofstream outfile;
    outfile.open(data_folder + "/cross_sections/total_" + projectile + "_" + target + "_" + sf_type + ".out");

    for (int ei = 0; ei < NE; ei++) {
        double E = pc->GeV * std::pow(10, logemin + ei * dE);
        double _xs;
        std::cout << "E [GeV] = " << E / pc->GeV << std::endl;
        if (E / pc->GeV > 0) {
              _xs = xs->TotalXS(E);
        }
        else {
            _xs = -99;
        }
        // std::cout << ", sigma/E =  " << _xs * 1e38 / (E / pc->GeV) << std::endl;
        outfile << E << "," << _xs << "\n";
    }
    outfile.close();

    return 0;
}