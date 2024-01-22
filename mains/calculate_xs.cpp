#include "configuration.h"
#include "physconst.h"
#include "structure_function.h"
#include "cross_section.h"
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
    const std::string xs_type = argv[4]; // Which SFs to use total, light, charm, ..

    std::cout << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "Config Path: " << config_path << std::endl;
    std::cout << "Projectile: " << projectile << std::endl;
    std::cout << "Target: " << target << std::endl;
    std::cout << "XS Type: " << xs_type << std::endl;
    std::cout << "=============================================" << std::endl << std::endl;

    Configuration config = Configuration(config_path);
    config.Populate();
    
    std::string data_folder = "../data/" + config.general.unique_name;

    config.Set_Target(target);
    config.Set_Projectile(projectile);
    config.Set_SF_Type(xs_type);

    PhysConst* pc = new PhysConst();
    CrossSection* xs = new CrossSection(config);

    // load the three structure function fit files
    string f1 = data_folder + "/F1_" + projectile + "_" + target + "_" + xs_type + ".fits";
    string f2 = data_folder + "/F2_" + projectile + "_" + target + "_" + xs_type + ".fits";
    string f3 = data_folder + "/F3_" + projectile + "_" + target + "_" + xs_type + ".fits";

    // TOTAL XS
    double logemin = 1;
    double logemax = 10;
    int NE = 200;
    double dE = (logemax - logemin) / (NE-1);

    int Ny = 100;
    double logymin = -6;
    double logymax = 0;
    double dy = (logymax - logymin) / (Ny-1);

    // Make the cross sections folder if it doesn't exist
    boost::filesystem::path out_folder = data_folder + "/cross_sections/";
    if (!boost::filesystem::exists(out_folder)) {
        boost::filesystem::create_directories(out_folder);
    }

    if (config.SF.mass_scheme != "parton") {
        xs->Load_Structure_Functions(f1, f2, f3);
    }

    std::ofstream outfile;
    outfile.open(data_folder + "/cross_sections/total_" + projectile + "_" + target + "_" + xs_type + ".out");

    for (int ei = 0; ei < NE; ei++) {
        double E = pc->GeV * std::pow(10, logemin + ei * dE);
        double _xs = std::log10(xs->TotalXS(E));
        outfile << E << "," << _xs << "\n";
    }
    outfile.close();

    // ds/dy
    std::ofstream dsdy_outfile;
    dsdy_outfile.open(data_folder + "/cross_sections/dsdy_" + projectile + "_" + target + "_" + xs_type + ".out");
    for (int ei = 0; ei < NE; ei++) {
        double E = pc->GeV * std::pow(10, logemin + ei * dE);
        for (int yi = 0; yi < Ny; yi++) { // loop over y
            double y = std::pow(10, logymin + yi * dy);
            // double y = ymin + yi * dy;
            // std::cout << E / 1e9 << " " << y << std::endl;
            double _dxs = std::log10(xs->ds_dy(E, y));
            // dsdy_outfile << (E / pc->GeV) << "," << y << "," << _dxs << std::endl;
            dsdy_outfile << _dxs << std::endl;
        }
        // dsdy_outfile << "\n";
    }
    dsdy_outfile.close();

    return 0;
}