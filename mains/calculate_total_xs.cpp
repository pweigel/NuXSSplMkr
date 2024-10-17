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
    if (argc != 11) {
        std::cout << "Not enough/too many inputs!" << std::endl;
        std::cout << "Usage: calculate_total_xs CONFIG CURRENT PROJECTILE TARGET TYPE MODE REPLICA EMIN EMAX LEPTON_FLAVOR" << std::endl;
        return 1;
    }

    const std::string config_path = argv[1]; // Path to .json file containing configuration info
    const std::string current = argv[2]; // CC or NC
    const std::string projectile = argv[3]; // neutrino or antineutrino
    const std::string target = argv[4]; // proton or neutron
    const std::string xs_type = argv[5]; // Which SFs to use total, light, charm, ..
    const int mode = std::stoi(argv[6]);
    const unsigned int replica = std::stoi(argv[7]);
    // double logemin = std::log10(1e1);
    // double logemax = std::log10(5e12);
    // if (argc > 8) {
    const double logemin = std::stod(argv[8]);
    // }
    // if (argc > 9) {
    const double logemax = std::stod(argv[9]);
    // }
    const std::string lepton_flavor = argv[10];

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

    // double logemin = std::log10(1e1);
    // double logemax = std::log10(5e12);

    int NE = 120;
    double dE = (logemax - logemin) / (NE-1);
    const vector<double> EnuTab{5e1, 1e2, 2e2, 5e2, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4,
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
    boost::filesystem::path out_folder = data_folder + "/cross_sections/"+lepton_flavor;
    if (!boost::filesystem::exists(out_folder)) {
        boost::filesystem::create_directories(out_folder);
    }
    
    config.Set_Current(current);
    config.Set_Projectile(projectile);
    config.Set_Target(target);
    config.Set_SF_Type(xs_type);

    if (lepton_flavor == "electron") {config.Set_Lepton_Mass(pc->electron_mass);}
    else if (lepton_flavor == "muon") {config.Set_Lepton_Mass(pc->muon_mass);}
    else if (lepton_flavor == "tau") {config.Set_Lepton_Mass(pc->tau_mass);}
    else { return 1; }

    if (current == "NC") {config.Set_Lepton_Mass(0.);}

    config.Set_Mode(mode);

    PhaseSpace ps(config);
    ps.Print();

    CrossSection* xs = new CrossSection(config, ps);
    std::string outfilename = data_folder + "/cross_sections/"+lepton_flavor+"/total_" + current + "_" + projectile + "_" + target + "_" + xs_type + "."+std::to_string(mode)+".out";
    if (config.XS.enable_radiative_corrections) {
        std::cout << "Radiative corrections enabled!" << std::endl;
        xs->Load_InterpGrid(data_folder + "/cross_sections"+lepton_flavor+"/dsdxdy_" + current + "_" + projectile + "_" + target + "_" + xs_type + "."+std::to_string(mode) + ".out");
        outfilename = data_folder + "/cross_sections/"+lepton_flavor+"total_" + current + "_" + projectile + "_" + target + "_" + xs_type + "."+std::to_string(mode) + ".rc";
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

    // std::ofstream outfile_table;
    // outfilename = data_folder + "/cross_sections/total_" + current + "_" + projectile + "_" + target + "_" + xs_type + "."+std::to_string(mode) + ".table";
    // if (config.XS.enable_radiative_corrections) {
    //     outfilename = data_folder + "/cross_sections/total_" + current + "_" + projectile + "_" + target + "_" + xs_type+ "."+std::to_string(mode)+ ".rctable";
    // }
    // outfile_table.open(outfilename);
    // for (auto const& Enu : EnuTab) {
    //     double E = pc->GeV * Enu;
    //     double _xs;
    //     std::cout << "E [GeV] = " << E / pc->GeV << std::endl;
    //     _xs = xs->TotalXS(E);
    //     outfile_table << E << "," << _xs << "\n";
    // }

    // outfile_table.close();

    return 0;
}