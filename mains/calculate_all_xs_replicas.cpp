#include "configuration.h"
#include "physconst.h"
#include "structure_function.h"
#include "cross_section.h"
#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>

using namespace nuxssplmkr;

int main(int argc, char* argv[]){
    if (argc < 3) {
        std::cout << "Not enough inputs!" << std::endl;
        return 1;
    } else if (argc > 3) {
        std::cout << "Too many inputs!" << std::endl;
        return 1;
    }

    const std::string config_path = argv[1]; // Path to .json file containing configuration info
    const std::string replica_string = argv[2];
    const int replica = stoi(replica_string);

    string targets[] = {"proton", "neutron"};
    string projectiles[] = {"neutrino", "antineutrino"};
    string sf_types[] = {"light", "charm"};

    std::cout << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "Config Path: " << config_path << std::endl;
    std::cout << "Making all cross sections!" << std::endl;
    std::cout << "REPLICA = " << replica_string << std::endl;
    std::cout << "=============================================" << std::endl << std::endl;

    double logemin = 2;
    double logemax = 9;
    int NE = 200;
    double dE = (logemax - logemin) / (NE-1);

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
                    // std::cout << "E = " << E / pc->GeV << " GeV" << std::endl;
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

            }
        }
    }

    return 0;
}