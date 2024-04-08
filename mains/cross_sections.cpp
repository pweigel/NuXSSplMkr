
#include "configuration.h"
#include "physconst.h"
#include "structure_function.h"
#include "cross_section.h"
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>

using namespace nuxssplmkr;
namespace po = boost::program_options;

int main(int argc, char* argv[]){
    std::string config_path;
    std::string projectile;
    std::string target;
    std::string xs_type;
    po::options_description desc("Options");
    desc.add_options()
        ("help", "Print help message")
        ("config_path", po::value<std::string>(&config_path), "Path to configuration file")
        ("projectile", po::value<std::string>(&projectile), "'neutrino' or 'antineutrino'")
        ("target", po::value<std::string>(&target), "'proton' or 'neutron'")
        ("xs_type", po::value<std::string>(&xs_type), "Which flavors to consider")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }
    
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

    // load the three structure function fit files
    string f1 = data_folder + "/F1_" + projectile + "_" + target + "_" + xs_type + ".fits";
    string f2 = data_folder + "/F2_" + projectile + "_" + target + "_" + xs_type + ".fits";
    string f3 = data_folder + "/F3_" + projectile + "_" + target + "_" + xs_type + ".fits";

    // TOTAL XS
    int NE = 100;
    int Ny = 100;
    int Nx = 100;

    double logemin = 2;
    double logemax = 12;
    double dE = (logemax - logemin) / (NE-1);

    double ymin = -1;
    double ymax = 0;
    double dy = (ymax - ymin) / (Ny-1);

    double xmin = -5;
    double xmax = 0;
    double dx = (xmax - xmin) / (Nx-1);

    std::ofstream outfile;
    outfile.open(data_folder + "/cross_sections/dsdxdy_" + projectile + "_" + target + "_" + xs_type + ".out");

    xs->Load_Structure_Functions(f1, f2, f3);
    xs->Set_Lepton_Mass(pc->muon_mass);
    
    double E = 31440.4 * pc->GeV;
    // double _xs = std::log10(xs->TotalXS(E));
    // std::cout << _xs << std::endl;
    for (int yi = 0; yi < Ny; yi++) { // loop over y
        double y = pow(10, ymin + yi * dy);
        for (int xi = 0; xi < Nx; xi++) { // loop over x
            double x = pow(10, xmin + xi * dx);
            // std::cout << "x=" << x << ", " << "y=" << y << std::endl;
            double _dsdxdy;
            if (!xs->PhaseSpaceIsGood(x, y, E)) {
                _dsdxdy = -99;
            } else {
                _dsdxdy = std::log10(xs->ds_dxdy(E, x, y));
            }
             
            outfile << _dsdxdy;

            if ( !((yi == Ny - 1) && (xi == Nx - 1)) ) {
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