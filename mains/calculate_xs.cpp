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
    if (argc < 7) {
        std::cout << "Not enough inputs!" << std::endl;
        return 1;
    } else if (argc > 7) {
        std::cout << "Too many inputs!" << std::endl;
        return 1;
    }
    const std::string config_path = argv[1]; // Path to .json file containing configuration info
    const std::string projectile = argv[2]; // neutrino or antineutrino
    const std::string target = argv[3]; // proton or neutron
    const std::string sf_type = argv[4]; // Which SFs to use total, light, charm, ..
    const int mode = std::stoi(argv[5]);
    const unsigned int replica = std::stoi(argv[6]);

    std::cout << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "Config Path: " << config_path << std::endl;
    std::cout << "Projectile: " << projectile << std::endl;
    std::cout << "Target: " << target << std::endl;
    std::cout << "SF Type: " << sf_type << std::endl;
    std::cout << "Mode: " << mode << std::endl;
    std::cout << "Replica: " << replica << std::endl;
    std::cout << "=============================================" << std::endl << std::endl;

    double logemin = 1;
    double logemax = 4;

    int NE = 110;
    double dE = (logemax - logemin) / (NE-1);
    const vector<double> EnuTab{1e2, 2e2, 5e2, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4,
      1e5, 2e5, 5e5, 1e6, 2e6, 5e6, 1e7, 2e7, 5e7, 1e8, 2e8,
      5e8, 1e9, 2e9, 5e9, 1e10, 2e10, 5e10, 1e11, 2e11, 5e11,
      1e12, 2e12, 5e12};

    PhysConst* pc = new PhysConst();

    Configuration config = Configuration(config_path);
    config.Populate();
    config.Set_Replica(replica);
    std::string data_folder = config.general.data_path + "/" + config.general.unique_name + "/replica_" + std::to_string(config.pdf_info.replica);
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
    config.Set_Mode(mode);
    std::cout << "Mode set to: " << mode << std::endl;

    PhaseSpace ps(config);
    ps.Print();

    CrossSection* xs = new CrossSection(config, ps);
    // if (config.SF.mass_scheme != "parton") {
    //     xs->Load_Structure_Functions(data_folder);
    // }
    xs->Set_Lepton_Mass(pc->muon_mass);

    // total
    std::ofstream outfile;
    outfile.open(data_folder + "/cross_sections/total_" + projectile + "_" + target + "_" + sf_type + ".out");

    // smooth E dist
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
        outfile << E
         << "," << _xs << "\n";
    }

    // std::ofstream outfile_table;
    // outfile_table.open(data_folder + "/cross_sections/total_" + projectile + "_" + target + "_" + sf_type + ".table");
    // for (auto const& Enu : EnuTab) {
    //     double E = pc->GeV * Enu;
    //     double _xs;
    //     std::cout << "E [GeV] = " << E / pc->GeV << std::endl;
    //     if (E / pc->GeV > 0) {
    //           _xs = xs->TotalXS(E);
    //     }
    //     else {
    //         _xs = -99;
    //     }
    //     outfile_table << E
    //      << "," << _xs << "\n";
    // }

    outfile.close();
    // outfile_table.close();

    return 0;
}