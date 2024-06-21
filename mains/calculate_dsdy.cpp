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
    if (argc != 7) {
        std::cout << "Not enough/too many inputs!" << std::endl;
        std::cout << "Usage: calculate_dsdy CONFIG PROJECTILE TARGET TYPE MODE REPLICA" << std::endl;
        return 1;
    }

    const std::string config_path = argv[1];
    const std::string projectile = argv[2]; // neutrino or antineutrino
    const std::string target = argv[3]; // proton or neutron
    const std::string xs_type = argv[4]; // Which SFs to use total, light, charm, ..
    const int mode = std::stoi(argv[5]);
    const unsigned int replica = std::stoi(argv[6]);

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

    config.Set_Projectile(projectile);
    config.Set_Target(target);
    config.Set_SF_Type(xs_type);
    config.Set_Lepton_Mass(pc->muon_mass);
    config.Set_Mode(mode);

    PhaseSpace ps(config);
    ps.Print();

    CrossSection* xs = new CrossSection(config, ps);

    int NE = 110;
    int Ny = 100;

    double logemin = std::log10(5e1);
    double logemax = std::log10(5e12);
    double dE = (logemax - logemin) / (NE-1);

    const vector<double> EnuTab{5e1, 1e2, 2e2, 5e2, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4,
      1e5, 2e5, 5e5, 1e6, 2e6, 5e6, 1e7, 2e7, 5e7, 1e8, 2e8,
      5e8, 1e9, 2e9, 5e9, 1e10, 2e10, 5e10, 1e11, 2e11, 5e11,
      1e12, 2e12, 5e12};

    // const vector<double> EnuTab{1e1, 2e1, 5e1, 1e2, 2e2, 5e2, 1e3, 2e3, 5e3, 1e4};

    // const vector<double> yvals{
    //     1e-6, 2e-6, 5e-6, 1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 
    //     1e-2, 2e-2, 5e-2, 8e-2, 1e-1, 1.2e-1, 1.5e-1, 1.8e-1, 2e-1, 2.2e-1, 2.5e-1, 2.8e-1,
    //     3e-1, 3.2e-1, 3.5e-1, 3.8e-1, 4e-1, 4.2e-1, 4.5e-1, 4.8e-1, 5e-1, 5.2e-1, 5.5e-1, 5.8e-1,
    //     6e-1, 6.2e-1, 6.5e-1, 6.8e-1, 7e-1, 7.2e-1, 7.5e-1, 7.8e-1, 8e-1, 8.2e-1, 8.5e-1, 8.8e-1,
    //     9e-1, 9.2e-1, 9.5e-1, 9.8e-1, 9.9e-1};

    double logymin = -4;
    double logymax = -1;
    double dy = (logymax - logymin) / Ny;

    int Nylin = 100;

    // ds/dy
    std::ofstream dsdy_outfile;
    dsdy_outfile.open(data_folder + "/cross_sections/dsdy_" + projectile + "_" + target + "_" + xs_type + ".out");

    // Get the energy, inelasticity values and put them in the header
    std::vector<double> energy_values;
    dsdy_outfile << "E";
    for (const auto E : EnuTab) { // loop over E
        dsdy_outfile << "," << E * pc->GeV;
    }
    dsdy_outfile << std::endl;

    std::vector<double> yvals;
    dsdy_outfile << "y";
    for (int i = 0; i < Ny; i++) {
        double y = pow(10, logymin + i * dy);
        dsdy_outfile << "," << y;
        yvals.push_back(y);
    }
    dy = (1.0 - pow(10, logymax))/Nylin;
    for (int i = 0; i < Nylin; i++) {
        double y = pow(10, logymax) + i*dy;
        dsdy_outfile << "," << y;
        yvals.push_back(y);
    }
    dsdy_outfile << std::endl;

    for (const auto E : EnuTab) { // loop over E
        std::cout << "E = " << E << " [GeV]" << std::endl;
        for (const auto y : yvals) {
            double _dxs = xs->ds_dy(E*pc->GeV, y);
            dsdy_outfile << _dxs << " ";
        }
        dsdy_outfile << std::endl;
    }
    dsdy_outfile.close();

    return 0;
}