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
    if (argc != 8) {
        std::cout << "Not enough/too many inputs!" << std::endl;
        std::cout << "Usage: calculate_total_xs CONFIG CURRENT PROJECTILE TARGET TYPE MODE REPLICA" << std::endl;
        return 1;
    }

    const std::string config_path = argv[1]; // Path to .json file containing configuration info
    const std::string current = argv[2]; // CC or NC
    const std::string projectile = argv[3]; // neutrino or antineutrino
    const std::string target = argv[4]; // proton or neutron
    const std::string xs_type = argv[5]; // Which SFs to use total, light, charm, ..
    const int mode = std::stoi(argv[6]);
    const unsigned int replica = std::stoi(argv[7]);

    std::cout << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "Config Path: " << config_path << std::endl;
    std::cout << "Current: " << current << std::endl;
    std::cout << "Projectile: " << projectile << std::endl;
    std::cout << "Target: " << target << std::endl;
    std::cout << "SF Type: " << xs_type << std::endl;
    std::cout << "Mode: " << mode << std::endl;
    std::cout << "Replica: " << replica << std::endl;
    std::cout << "=============================================" << std::endl << std::endl;

    int NE = 551;
    int Nz = 551;

    double logemin = std::log10(1e1);
    double logemax = std::log10(1e12);
    double dE = (logemax - logemin) / (NE-1);

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
    
    config.Set_Current(current);
    config.Set_Projectile(projectile);
    config.Set_Target(target);
    config.Set_SF_Type(xs_type);
    config.Set_Lepton_Mass(pc->muon_mass);
    config.Set_Mode(mode);

    PhaseSpace ps(config);
    ps.Print();

    CrossSection* xs = new CrossSection(config, ps);
    std::string outfilename = data_folder + "/cross_sections/nusquids_sigma_" + current + "_" + projectile + "_" + target + "_" + xs_type + "."+std::to_string(mode)+".out";
    if (config.XS.enable_radiative_corrections) {
        std::cout << "Radiative corrections enabled!" << std::endl;
        xs->Load_InterpGrid(data_folder + "/cross_sections/dsdxdy_" + current + "_" + projectile + "_" + target + "_" + xs_type + "."+std::to_string(mode) + ".out");
        outfilename = data_folder + "/cross_sections/nusquids_sigma_" + current + "_" + projectile + "_" + target + "_" + xs_type + "."+std::to_string(mode) + ".rc";
    }

    std::ofstream outfile;
    outfile.open(outfilename);

    // smooth E dist
    for (int ei = 0; ei < NE; ei++) {
        double E = pc->GeV * std::pow(10, logemin + ei * dE);
        double _xs;
        std::cout << "E [GeV] = " << E / pc->GeV << std::endl;
        _xs = xs->TotalXS(E);
        outfile << E << "," << _xs << "\n";
    }
    outfile.close();

    return 0;
}