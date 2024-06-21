
#include "NuXSSplMkr/configuration.h"
#include "NuXSSplMkr/physconst.h"
#include "NuXSSplMkr/phase_space.h"
#include "NuXSSplMkr/structure_function.h"
#include "NuXSSplMkr/cross_section.h"
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>

using namespace nuxssplmkr;
namespace po = boost::program_options;

int main(int argc, char* argv[]){

   if (argc < 6) {
        std::cout << "Not enough inputs." << std::endl;
        return 0;
    }
    const std::string config_path = argv[1];
    const std::string projectile = argv[2]; // neutrino or antineutrino
    const std::string target = argv[3]; // proton or neutron
    const std::string xs_type = argv[4]; // Which SFs to use total, light, charm, ..
    const int mode = std::stoi(argv[5]);
    const unsigned int replica = std::stoi(argv[6]);

    std::cout << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "Config Path: " << config_path << std::endl;
    std::cout << "Projectile: " << projectile << std::endl;
    std::cout << "Target: " << target << std::endl;
    std::cout << "XS Type: " << xs_type << std::endl;
    std::cout << "Replica: " << replica << std::endl;
    std::cout << "Mode: " << mode << std::endl;
    std::cout << "=============================================" << std::endl << std::endl;

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

    config.Set_Projectile(projectile);
    config.Set_Target(target);
    config.Set_SF_Type(xs_type);
    config.Set_Mode(mode);

    PhysConst* pc = new PhysConst();
    PhaseSpace ps(config);
    ps.Print();

    CrossSection* xs = new CrossSection(config, ps);

    // int NE = 110;
    int Ny = 300;
    int Nx = 300;

    double logemin = 1;
    double logemax = 12;
    int NE = 1;
    // double dE = (logemax - logemin) / (NE-1);
    double dE = 1.;

    double logymin = -3;
    double logymax = 0;
    double dy = (logymax - logymin) / (Ny-1);

    double logxmin = -3;
    double logxmax = 0;
    double dx = (logxmax - logxmin) / (Nx-1);

    std::ofstream outfile;
    outfile.open(data_folder + "/cross_sections/dsdxdy_" + projectile + "_" + target + "_" + xs_type + ".out");

    // xs->Load_Structure_Functions(f1, f2, f3);
    // xs->Load_Structure_Functions(data_folder);
    xs->Set_Lepton_Mass(pc->muon_mass);
    
    // Get the energy, inelasticity values and put them in the header
    std::vector<double> energy_values;
    outfile << "E";
    for (int ei = 0; ei < NE; ei++) {
        double E = pc->GeV * std::pow(10, logemin + ei * dE);
        energy_values.push_back(E);
        outfile << "," << E;
    }
    outfile << std::endl;

    std::vector<double> inelasticity_values;
    outfile << "y";
    for (int yi = 0; yi < Ny; yi++) {
        double y = std::pow(10, logymin + yi * dy);
        inelasticity_values.push_back(y);
        outfile << "," << y;
    }
    outfile << std::endl;

    std::vector<double> x_values;
    outfile << "x";
    for (int xi = 0; xi < Nx; xi++) {
        double x = std::pow(10, logxmin + xi * dx);
        x_values.push_back(x);
        outfile << "," << x;
    }
    outfile << std::endl;

    for (int ei = 0; ei < NE; ei++) {
        double E = energy_values[ei];
        xs->Set_Neutrino_Energy(E);
        for (int yi = 0; yi < Ny; yi++) { // loop over y
            double y = inelasticity_values[yi];
            for (int xi = 0; xi < Nx; xi++) { // loop over x
                double x = x_values[xi];
                double _dsdxdy;
                bool ps_valid = ps.Validate(E, x, y);

                if (!ps_valid) {
                    _dsdxdy = 0.0;
                } else {
                    _dsdxdy = xs->ds_dxdy(x, y);
                }
                
                outfile << _dsdxdy;
                if ( !(xi == Nx - 1) ) {
                    outfile << ",";
                }
            }
            outfile << std::endl;
        }
    }

    outfile.close();

    return 0;
}