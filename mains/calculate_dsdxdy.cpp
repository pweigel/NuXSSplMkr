
#include "configuration.h"
#include "physconst.h"
#include "phase_space.h"
#include "structure_function.h"
#include "cross_section.h"
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>

using namespace nuxssplmkr;
namespace po = boost::program_options;

int main(int argc, char* argv[]){

   if (argc < 2) {
        std::cout << "Not enough inputs." << std::endl;
        return 0;
    }
    const std::string config_path = argv[1];
    const std::string projectile = argv[2]; // neutrino or antineutrino
    const std::string target = argv[3]; // proton or neutron
    const std::string xs_type = argv[4]; // Which SFs to use total, light, charm, ..

    const int replica = 0; // TODO: make this an input

    // Create a new config w/ the filename
    std::cout << config_path << std::endl;
    Configuration config = Configuration(config_path);
    config.Populate();
    config.Set_Replica(replica);
    std::string data_folder = "../data/" + config.general.unique_name + "/replica_" + std::to_string(config.pdf.replica);
    std::cout << "Loading/saving data to: " << data_folder << std::endl;
    
    // Make the cross sections folder if it doesn't exist
    boost::filesystem::path out_folder = data_folder + "/cross_sections/";
    if (!boost::filesystem::exists(out_folder)) {
        boost::filesystem::create_directories(out_folder);
    }

    config.Set_Projectile(projectile);
    config.Set_Target(target);
    config.Set_SF_Type(xs_type);

    PhysConst* pc = new PhysConst();
    PhaseSpace ps(config);
    ps.Print();

    CrossSection* xs = new CrossSection(config, ps);

    // load the three structure function fit files
    string f1 = data_folder + "/F1_" + projectile + "_" + target + "_" + xs_type + ".fits";
    string f2 = data_folder + "/F2_" + projectile + "_" + target + "_" + xs_type + ".fits";
    string f3 = data_folder + "/F3_" + projectile + "_" + target + "_" + xs_type + ".fits";

    // TOTAL XS
    int NE = 200;
    int Ny = 100;
    int Nx = 100;

    double logemin = 2;
    double logemax = 9;
    double dE = (logemax - logemin) / (NE-1);

    double logymin = -7;
    double logymax = 0;
    double dy = (logymax - logymin) / (Ny-1);

    double logxmin = -7;
    double logxmax = 0;
    double dx = (logxmax - logxmin) / (Nx-1);

    std::ofstream outfile;
    outfile.open(data_folder + "/cross_sections/dsdxdy_" + projectile + "_" + target + "_" + xs_type + ".out");

    xs->Load_Structure_Functions(f1, f2, f3);
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
        for (int yi = 0; yi < Ny; yi++) { // loop over y
            double y = inelasticity_values[yi];
            for (int xi = 0; xi < Nx; xi++) { // loop over x
                double x = x_values[xi];
                double _dsdxdy;

                bool valid = ps.Validate(E, x, y);
                if (!valid) {
                    _dsdxdy = 0.0;
                } else {
                    _dsdxdy = std::log10(xs->ds_dxdy(E, x, y));
                }

                outfile << _dsdxdy;
                if ( !((yi == Ny - 1) && (xi == Nx - 1)) ) {
                    outfile << ",";
                }
            }
            outfile << std::endl;
        }
        outfile << std::endl;
    }

    outfile.close();

    return 0;
}