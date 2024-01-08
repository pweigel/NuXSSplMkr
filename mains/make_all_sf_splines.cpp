#include "configuration.h"
#include "physconst.h"
#include "structure_function.h"
#include <boost/filesystem.hpp>

int main(int argc, char* argv[]){
    if (argc < 2) {
        std::cout << "Not enough inputs!" << std::endl;
        return 1;
    } else if (argc > 2) {
        std::cout << "Too many inputs!" << std::endl;
        return 1;
    }

    const std::string config_path = argv[1];

    string targets[] = {"proton", "neutron"};
    string projectiles[] = {"neutrino", "antineutrino"};
    string sf_types[] = {"total", "charm"};

    std::cout << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "Config Path: " << config_path << std::endl;
    std::cout << "Making all splines!" << std::endl;
    std::cout << "=============================================" << std::endl << std::endl;

    // Create a new config w/ the filename
    nuxssplmkr::Configuration config = nuxssplmkr::Configuration(config_path);
    config.Populate();  // Populate the SF info

    for (const string &sf_type : sf_types) {
        config.Set_SF_Type(sf_type);
        std::cout << "SFType set to: " << sf_type << std::endl;

        for (const string &projectile : projectiles) {
            config.Set_Projectile(projectile);
            std::cout << "Projectile set to: " << projectile << std::endl;

            for (const string &target : targets) {
                config.Set_Target(target);
                std::cout << "Target set to: " << target << std::endl;

                nuxssplmkr::StructureFunction sf = nuxssplmkr::StructureFunction(config);
                nuxssplmkr::PhysConst* pc = new nuxssplmkr::PhysConst();

                sf.Set_Lepton_Mass(pc->muon_mass);
                
                // std::string _out_folder = "../data/" + config.sf_info.pdfset + "_" + config.sf_info.mass_scheme + "_pto" + to_string(config.sf_info.perturbative_order);
                std::string _out_folder = "../data/" + config.sf_info.pdfset + "_" + config.sf_info.mass_scheme + "_pto" + to_string(config.sf_info.perturbative_order) + "_" + config.sf_info.small_x_order;
                boost::filesystem::path out_folder = _out_folder;
                if (!boost::filesystem::exists(out_folder)) {  // make the out_folder if it does not exist
                    boost::filesystem::create_directories(out_folder);
                }
                
                sf.InitializeAPFEL();
                sf.BuildSplines(out_folder.string()); // Photospline
            }
        }

    }

    return 0;
}