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
    if (argc != 8) {
        std::cout << "Not enough/too many inputs!" << std::endl;
        std::cout << "Usage: calculate_dsdy CONFIG CURRENT PROJECTILE TARGET TYPE MODE REPLICA" << std::endl;
        return 1;
    }

    const std::string config_path = argv[1];
    const std::string current = argv[2]; // CC or NC
    const std::string projectile = argv[3]; // neutrino or antineutrino
    const std::string target = argv[4]; // proton or neutron
    const std::string xs_type = argv[5]; // Which SFs to use total, light, charm, ..
    const int mode = std::stoi(argv[6]);
    const unsigned int replica = std::stoi(argv[7]);
    // const bool use_rc = !!(std::stoi(argv[7]));

    // Create a new config w/ the filename
    std::cout << config_path << std::endl;
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
    
    config.Set_Current(current);
    config.Set_Projectile(projectile);
    config.Set_Target(target);
    config.Set_SF_Type(xs_type);
    config.Set_Lepton_Mass(pc->muon_mass);
    config.Set_Mode(mode);

    PhaseSpace ps(config);
    ps.Print();

    CrossSection* xs = new CrossSection(config, ps);
    std::string outfilename = data_folder + "/cross_sections/nusquids_dsdy_" + current + "_" + projectile + "_" + target + "_" + xs_type + "."+std::to_string(mode) + ".out";
    if (config.XS.enable_radiative_corrections) {
        std::cout << "Radiative corrections enabled!" << std::endl;
        xs->Load_InterpGrid(data_folder + "/cross_sections/dsdxdy_" + current + "_" + projectile + "_" + target + "_" + xs_type + "."+std::to_string(mode) + ".out");
        outfilename = data_folder + "/cross_sections/nusquids_dsdy_" + current + "_" + projectile + "_" + target + "_" + xs_type + "."+std::to_string(mode) + ".rc";
    }

    int NE = 551;
    int Nz = 551;

    // y = (Ein - Eout) / Ein
    // Eout = Ein * (1-y)
    // z = (Eout - Eout,min)/(Ein - Eout,min)
    // z = (Ein(1-y) - Eout,min)/(Ein - Eout,min)
    // z (Ein - Eout,min) = Ein(1-y) - Eout,min
    // z (Ein - Eout,min) + Eout,min = Ein(1-y)
    // (z (Ein - Eout,min) + Eout,min) / Ein = 1-y
    // y = 1 - (z (Ein - Eout,min) + Eout,min) / Ein
    // Alfonso:
    // y = (1-z)(ei-10)/ei

    double logemin = std::log10(1e1);
    double logemax = std::log10(1e12);
    double dE = (logemax - logemin) / (NE-1);

    double emin = pow(10.0, logemin);
    double emax = pow(10.0, logemax);

    double zmin = 0.0;
    double zmax = 1.0;
    double dz = (zmax - zmin) / (Nz - 1);

    // ds/dy
    std::ofstream dsdy_outfile;
    dsdy_outfile.open(outfilename);

    // Get the energy, z values and put them in the header
    std::vector<double> energy_values;
    dsdy_outfile << "E";
    for (int i = 0; i < NE; i++) { // loop over E
        double logE = logemin + i * dE;
        double E = pow(10.0, logE);
        energy_values.push_back(E);
        dsdy_outfile << "," << E * pc->GeV;
    }
    dsdy_outfile << std::endl;

    std::vector<double> zvals;
    dsdy_outfile << "z";
    for (int i = 0; i < Nz; i++) {
        double z = zmin + i * dz;
        dsdy_outfile << "," << z;
        zvals.push_back(z);
    }
    dsdy_outfile << std::endl;

    for (const auto E : energy_values) { // loop over E
        std::cout << "E = " << E << " [GeV]" << std::endl;
        for (int i = 0; i < Nz; i++) {
            double z = zvals[i];
            double y = (1 - z) * (E - emin) / E;
            // use alfonso's modification
            if (y == 1.) y -= 1e-4;
            if (y == 0.) {
                y = 5.0 / E;
                double z_prev = zvals[i-1]; 
                double y_prev = (1 - z_prev) * (E - emin) / E;
                if (y > y_prev) y = 1e-4;
            }
            // std::cout << y << std::endl;
            double _dxs = xs->ds_dy(E*pc->GeV, y);
            if (_dxs == 0.0) {
                _dxs = 1e-50;
            }
            dsdy_outfile << std::log10(_dxs) << " ";
        }
        dsdy_outfile << std::endl;
    }
    dsdy_outfile.close();

    return 0;
}