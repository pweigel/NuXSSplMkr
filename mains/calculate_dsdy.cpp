#include "configuration.h"
#include "physconst.h"
#include "structure_function.h"
#include "cross_section.h"
#include "tools.h"
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
    const std::string projectile = argv[2]; // neutrino or antineutrino
    const std::string target = argv[3]; // proton or neutron
    const std::string xs_type = argv[4]; // Which SFs to use total, light, charm, ..

    // Create a new config w/ the filename
    std::cout << config_path << std::endl;
    Configuration config = Configuration(config_path);
    config.Populate();  // Populate the SF info
    
    std::string data_folder = "../data/" + config.general.unique_name;

    config.Set_Projectile(projectile);
    config.Set_Target(target);
    config.Set_SF_Type(xs_type);

    PhysConst* pc = new PhysConst();
    CrossSection* xs = new CrossSection(config);

    // load the three structure function fit files
    string f1 = data_folder + "/F1_" + projectile + "_" + target + "_" + xs_type + ".fits";
    string f2 = data_folder + "/F2_" + projectile + "_" + target + "_" + xs_type + ".fits";
    string f3 = data_folder + "/F3_" + projectile + "_" + target + "_" + xs_type + ".fits";

    std::vector<string> fns = {f1, f2, f3};
    for (string fn : fns) {
        if (!fexists(fn)) {
            std::cout << "Could not find file: " << fn << std::endl;
        } else {
            std::cout << "Found file: " << fn << std::endl;
        }
    }

    int NE = 200;
    int Ny = 100;

    double logemin = 1;
    double logemax = 11;
    double dE = (logemax - logemin) / (NE-1);

    double logymin = -6;
    double logymax = 0;
    double dy = (logymax - logymin) / Ny;
    xs->Set_Lepton_Mass(pc->muon_mass);
    if (config.SF.mass_scheme != "parton") {
        xs->Load_Structure_Functions(f1, f2, f3);
    }

    // ds/dy
    std::ofstream dsdy_outfile;
    dsdy_outfile.open(data_folder + "/cross_sections/dsdy_" + projectile + "_" + target + "_" + xs_type + ".out");
    for (int ei = 0; ei < NE; ei++) {
        double E = pc->GeV * std::pow(10, logemin + ei * dE);
        std::cout << "E = " << E / (pc->GeV) << " GeV" << std::endl;
        for (int yi = 0; yi < Ny; yi++) { // loop over y
            double y = std::pow(10, logymin + yi * dy);
            double _dxs = std::log10(xs->ds_dy(E, y));
            dsdy_outfile << _dxs << std::endl;
        }
        // dsdy_outfile << "\n";
    }
    dsdy_outfile.close();

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