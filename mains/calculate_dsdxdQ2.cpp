
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
    const std::string projectile = argv[2]; // neutrino or antineutrino
    const std::string target = argv[3]; // proton or neutron
    const std::string xs_type = argv[4]; // Which SFs to use total, light, charm, ..

    // Create a new config w/ the filename
    std::cout << config_path << std::endl;
    Configuration config = Configuration(config_path);
    config.Populate();  // Populate the SF info
    
    std::string data_folder = "../data/" + config.general.unique_name;
    // Make the cross sections folder if it doesn't exist
    boost::filesystem::path out_folder = data_folder + "/cross_sections/";
    if (!boost::filesystem::exists(out_folder)) {
        boost::filesystem::create_directories(out_folder);
    }

    config.Set_Projectile(projectile);
    config.Set_Target(target);
    config.Set_SF_Type(xs_type);

    PhysConst* pc = new PhysConst();
    CrossSection* xs = new CrossSection(config);

    std::cout << "TEST" << std::endl;

    // load the three structure function fit files
    string f1 = data_folder + "/F1_" + projectile + "_" + target + "_" + xs_type + ".fits";
    string f2 = data_folder + "/F2_" + projectile + "_" + target + "_" + xs_type + ".fits";
    string f3 = data_folder + "/F3_" + projectile + "_" + target + "_" + xs_type + ".fits";

    // TOTAL XS
    int NQ2 = 100;
    int Ny = 100;
    int Nx = 100;

    double Q2min = 0;
    double Q2max = 4;
    double dQ2 = (Q2max - Q2min) / (NQ2-1);

    double xmin = -4;
    double xmax = 0;
    double dx = (xmax - xmin) / (Nx-1);

    std::ofstream outfile;
    outfile.open(data_folder + "/cross_sections/dsdxdQ2_" + projectile + "_" + target + "_" + xs_type + ".out");

    xs->Load_Structure_Functions(f1, f2, f3);
    xs->Set_Lepton_Mass(pc->muon_mass);
    
    double E = 1000.0 * pc->GeV;
    double _xs = std::log10(xs->TotalXS(1.0 * pc->PeV));
    double s = 2 * 0.938 * E;

    std::cout << _xs << std::endl;
    for (int qi = 0; qi < NQ2; qi++) { // loop over Q2
        double Q2 = pow(10, Q2min + qi * dQ2) * SQ(pc->GeV);
        for (int xi = 0; xi < Nx; xi++) { // loop over x
            double x = pow(10, xmin + xi * dx);
            xs->Set_Neutrino_Energy(E);
            // std::cout << "x=" << x << ", " << "y=" << y << std::endl;
            double _dsdxdQ2;

            // a test:
            if (Q2 > s && xs->PhaseSpaceIsGood_Q2(x, Q2, E)) {
                std::cout << Q2 / SQ(pc->GeV) << ", " << x << std::endl;
            }

            if (!xs->PhaseSpaceIsGood_Q2(x, Q2, E)) {
                _dsdxdQ2 = -99;
            } else {
                _dsdxdQ2 = std::log10(xs->ds_dxdQ2(E, x, Q2));
            }

            std::cout << Q2 / SQ(pc->GeV) << ", " << x << ", " << _dsdxdQ2 << std::endl;
             
            outfile << _dsdxdQ2;

            if ( !((qi == NQ2 - 1) && (xi == Nx - 1)) ) {
                outfile << ",";
            }
        }
    }

    outfile.close();

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