
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
    if (argc != 7) {
        std::cout << "Not enough/too many inputs!" << std::endl;
        std::cout << "Usage: calculate_dsdxdy CONFIG PROJECTILE TARGET TYPE MODE REPLICA" << std::endl;
        return 1;
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

    PhysConst* pc = new PhysConst();

    config.Set_Projectile(projectile);
    config.Set_Target(target);
    config.Set_SF_Type(xs_type);
    config.Set_Lepton_Mass(pc->muon_mass);
    config.Set_Mode(mode);

    PhaseSpace ps(config);
    ps.Print();

    CrossSection* xs = new CrossSection(config, ps);
    std::string outfilename = data_folder + "/cross_sections/dsdxdy_" + projectile + "_" + target + "_" + xs_type + ".out";
    if (config.XS.enable_radiative_corrections) {
        std::cout << "Radiative corrections enabled!" << std::endl;
        xs->Load_InterpGrid(data_folder + "/cross_sections/dsdxdy_" + projectile + "_" + target + "_" + xs_type + ".out");
        outfilename = data_folder + "/cross_sections/dsdxdy_" + projectile + "_" + target + "_" + xs_type + ".rc2";
    }
    
    int NE = 130;
    int Ny = 255;
    int Nylin = 100;
    int Nx = 200;
    int Nxlin = 100;

    double logemin = 1;
    double logemax = 13;
    double dE = (logemax - logemin) / (NE-1);

    double logymin = -16;
    double logymax = -1;
    double dy = (logymax - logymin) / (Ny);

    double logxmin = -12;
    double logxmax = -1;
    double dx = (logxmax - logxmin) / (Nx);

    std::ofstream outfile;
    outfile.open(outfilename);

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
    for (int i = 0; i < Ny; i++) {
        double y = pow(10, logymin + i * dy);
        outfile << "," << y;
        inelasticity_values.push_back(y);
    }
    dy = (1.0 - pow(10, logymax))/Nylin;
    for (int i = 0; i < Nylin; i++) {
        double y = pow(10, logymax) + i*dy;
        outfile << "," << y;
        inelasticity_values.push_back(y);
    }
    outfile << std::endl;

    std::vector<double> x_values;
    outfile << "x";
    for (int xi = 0; xi < Nx; xi++) {
        double x = std::pow(10, logxmin + xi * dx);
        x_values.push_back(x);
        outfile << "," << x;
    }
    dx = (1.0 - pow(10, logxmax))/Nxlin;
    for (int i = 0; i < Nylin; i++) {
        double x = pow(10, logxmax) + i*dx;
        x_values.push_back(x);
        outfile << "," << x;
    }
    outfile << std::endl;

    for (int ei = 0; ei < NE; ei++) {
        double E = energy_values[ei];
        xs->Set_Neutrino_Energy(E);
        std::cout << E/1e9 << std::endl;
        for (int yi = 0; yi < Ny+Nylin; yi++) { // loop over y
            double y = inelasticity_values[yi];
            for (int xi = 0; xi < Nx+Nxlin; xi++) { // loop over x
                double x = x_values[xi];
                double _dsdxdy;
                bool ps_valid = ps.Validate(E, x, y);

                if (!ps_valid) {
                    _dsdxdy = 0.0;
                } else {
                    _dsdxdy = xs->ds_dxdy(x, y);
                    double rc = 0.0;
                    if (config.XS.enable_radiative_corrections) {
                        rc = xs->rc_dsdxdy(E, x, y, _dsdxdy);
                        // rc = xs->rc_bardin(E, x, y);
                        // std::cout << rc << std::endl;
                    }
                    // if (rc / _dsdxdy < -1) {
                    //     double zmin = 1 - y + x * y;
                    //     double xhat = x*y / (zmin+y-1);
                    //     double yhat = (zmin+y-1)/zmin;
                    //     std::cout << scientific << setprecision(12);
                    //     std::cout << E/1e9 << ", " << x << ", " << y << " (" << zmin << ", " << xhat << ", " << yhat << ")" << ": " << _dsdxdy << ", " << rc << std::endl;
                    // }
                    _dsdxdy = _dsdxdy + rc;
                    // _dsdxdy *= (1.0 + rc);

                }
                
                outfile << _dsdxdy;
                if ( !(xi == Nx + Nxlin - 1) ) {
                    outfile << ",";
                }
            }
            outfile << std::endl;
        }
    }

    outfile.close();

    return 0;
}