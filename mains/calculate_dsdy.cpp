#include "NuXSSplMkr/configuration.h"
#include "NuXSSplMkr/physconst.h"
#include "NuXSSplMkr/phase_space.h"
#include "NuXSSplMkr/structure_function.h"
#include "NuXSSplMkr/cross_section.h"
#include "NuXSSplMkr/tools.h"
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

    int NE = 200;
    int Ny = 100;

    double logemin = 1;
    double logemax = 9;
    double dE = (logemax - logemin) / (NE-1);

    double logymin = -6;
    double logymax = 0;
    double dy = (logymax - logymin) / Ny;

    xs->Set_Lepton_Mass(pc->muon_mass);
    xs->Load_Structure_Functions(f1, f2, f3);

    // ds/dy
    std::ofstream dsdy_outfile;
    dsdy_outfile.open(data_folder + "/cross_sections/dsdy_" + projectile + "_" + target + "_" + xs_type + ".out");

    // Get the energy, inelasticity values and put them in the header
    std::vector<double> energy_values;
    dsdy_outfile << "E";
    for (int ei = 0; ei < NE; ei++) {
        double E = pc->GeV * std::pow(10, logemin + ei * dE);
        energy_values.push_back(E);
        dsdy_outfile << "," << E;
    }
    dsdy_outfile << std::endl;

    std::vector<double> inelasticity_values;
    dsdy_outfile << "y";
    for (int yi = 0; yi < Ny; yi++) {
        double y = std::pow(10, logymin + yi * dy);
        inelasticity_values.push_back(y);
        dsdy_outfile << "," << y;
    }
    dsdy_outfile << std::endl;

    for (int ei = 0; ei < NE; ei++) { // loop over E
        double E = energy_values[ei];
        std::cout << "E = " << E / (pc->GeV) << " GeV" << std::endl;
        for (int yi = 0; yi < Ny; yi++) { // loop over y
            double y = inelasticity_values[yi];
            double _dxs = xs->ds_dy(E, y);

            dsdy_outfile << _dxs;
            if ( !(yi == Ny - 1) ) {
                dsdy_outfile << ",";
            }
        }
        dsdy_outfile << std::endl;
    }
    dsdy_outfile.close();



    return 0;
}