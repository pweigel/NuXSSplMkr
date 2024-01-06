#include "configuration.h"
#include "physconst.h"
#include "structure_function.h"
#include "cross_section.h"
#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>

using namespace nuxssplmkr;

int main(int argc, char* argv[]){
    if (argc < 2) {
        std::cout << "Not enough inputs." << std::endl;
        return 0;
    }
    const std::string config_path = argv[1];
    const std::string target = argv[2]; // proton or neutron
    const std::string projectile = argv[3]; // neutrino or antineutrino
    const std::string xs_type = argv[4]; // Which SFs to use total, light, charm, ..

    // Create a new config w/ the filename
    std::cout << config_path << std::endl;
    Configuration config = Configuration(config_path);
    config.Populate();  // Populate the SF info
    
    std::string data_folder = "../data/" + config.sf_info.pdfset + "_" + config.sf_info.mass_scheme + "_pto" + to_string(config.sf_info.perturbative_order);

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
    int NE = 100;
    int Ny = 1000;

    double logemin = 2;
    double logemax = 9;
    double dE = (logemax - logemin) / (NE-1);

    double ymin = 1e-4;
    double ymax = 1;
    double dy = (ymax - ymin) / Ny;

    std::ofstream outfile;
    outfile.open(data_folder + "/cross_sections/dsdy_" + target + "_" + xs_type + ".out");
    xs->Load_Structure_Functions(f1, f2, f3);
    double E = 1000.0 * pc->GeV;
    for (int yi = 0; yi < Ny; yi++) { // loop over y
        double y = ymin + yi * dy;
        outfile << y << "," << std::log10(xs->ds_dy(E, y)) << std::endl;
    }

    // SINGLE DIFFERENTIAL XS, linear y
    // for (int ei = 0; ei < NE; ei++) { // loop over E
    //     double E = pc->GeV * std::pow(10, logemin + ei * dE);
    //     for (int yi = 0; yi < Ny; yi++) { // loop over y
    //         double y = ymin + yi * dy;
    //         outfile << std::log10(xs->ds_dy(E, y));
    //         if (yi != Ny - 1) {
    //             outfile << ",";
    //         }
    //     }
    //     outfile << "\n";
    // }

    return 0;
}