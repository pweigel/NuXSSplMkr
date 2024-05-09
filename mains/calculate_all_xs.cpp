#include "NuXSSplMkr/configuration.h"
#include "NuXSSplMkr/physconst.h"
#include "NuXSSplMkr/structure_function.h"
#include "NuXSSplMkr/cross_section.h"
#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>

using namespace nuxssplmkr;

int main(int argc, char* argv[]){
    if (argc < 2) {
        std::cout << "Not enough inputs!" << std::endl;
        return 1;
    } else if (argc > 2) {
        std::cout << "Too many inputs!" << std::endl;
        return 1;
    }
    const std::string config_path = argv[1]; // Path to .json file containing configuration info

    string targets[] = {"proton", "neutron"};
    string projectiles[] = {"neutrino", "antineutrino"};
    // string sf_types[] = {"total", "charm", "top"};
    // string targets[] = {"proton"};
    // string projectiles[] = {"neutrino"};
    // string sf_types[] = {"total"};
    string sf_types[] = {"light"};
    // string targets[] = {"proton", "neutron"};
    // string projectiles[] = {"neutrino"};
    // string sf_types[] = {"charm"};
    std::cout << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "Config Path: " << config_path << std::endl;
    std::cout << "Making all cross sections!" << std::endl;
    std::cout << "=============================================" << std::endl << std::endl;

    double logemin = 2;
    double logemax = 9;
    int NE = 200;
    double dE = (logemax - logemin) / (NE-1);

    // int Ny = 100;
    // double ymin = 1e-6;
    // double ymax = 1.0;
    // double dy = (ymax - ymin) / (Ny-1);
    // double logymin = -6;
    // double logymax = 0;
    // double dy = (logymax - logymin) / (Ny-1);

    PhysConst* pc = new PhysConst();

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
    
    for (const string &sf_type : sf_types) {
        config.Set_SF_Type(sf_type);
        std::cout << "SFType set to: " << sf_type << std::endl;
        for (const string &projectile : projectiles) {
            config.Set_Projectile(projectile);
            std::cout << "Projectile set to: " << projectile << std::endl;
            for (const string &target : targets) {
                config.Set_Target(target);
                std::cout << "Target set to: " << target << std::endl;

                string f1 = data_folder + "/F1_" + projectile + "_" + target + "_" + sf_type + ".fits";
                string f2 = data_folder + "/F2_" + projectile + "_" + target + "_" + sf_type + ".fits";
                string f3 = data_folder + "/F3_" + projectile + "_" + target + "_" + sf_type + ".fits";

                CrossSection* xs = new CrossSection(config);
                if (config.SF.mass_scheme != "parton") {
                    xs->Load_Structure_Functions(f1, f2, f3);
                }
                xs->Set_Lepton_Mass(pc->muon_mass);

                // total
                std::ofstream outfile;
                outfile.open(data_folder + "/cross_sections/total_" + projectile + "_" + target + "_" + sf_type + ".out");

                for (int ei = 0; ei < NE; ei++) {
                    double E = pc->GeV * std::pow(10, logemin + ei * dE);
                    std::cout << "E = " << E / pc->GeV << " GeV" << std::endl;
                    double _xs;
                    if (E / pc->GeV > 0) {
                         _xs = std::log10(xs->TotalXS(E));
                    }
                    else {
                        _xs = -99;
                    }
                    outfile << E << "," << _xs << "\n";
                }
                outfile.close();

                // ds/dy
                // std::ofstream dsdy_outfile;
                // dsdy_outfile.open(data_folder + "/cross_sections/dsdy_" + projectile + "_" + target + "_" + sf_type + ".out");
                // for (int ei = 0; ei < NE; ei++) {
                //     double E = pc->GeV * std::pow(10, logemin + ei * dE);
                //     std::cout << "E = " << E / pc->GeV << " GeV" << std::endl;
                //     for (int yi = 0; yi < Ny; yi++) { // loop over y
                //         double y = std::pow(10, logymin + yi * dy);
                //         // double y = ymin + yi * dy;
                //         // std::cout << E / 1e9 << " " << y << std::endl;
                //         double _dxs = std::log10(xs->ds_dy(E, y));
                //         // dsdy_outfile << (E / pc->GeV) << "," << y << "," << _dxs << std::endl;
                //         dsdy_outfile << _dxs << std::endl;
                //     }
                //     // dsdy_outfile << "\n";
                // }
                // dsdy_outfile.close();
            }
        }
    }



    return 0;
}