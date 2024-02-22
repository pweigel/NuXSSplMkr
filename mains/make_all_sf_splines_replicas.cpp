#include "configuration.h"
#include "physconst.h"
#include "structure_function.h"
#include <boost/filesystem.hpp>

int main(int argc, char* argv[]){
    if (argc < 2) {
        std::cout << "Not enough inputs!" << std::endl;
        return 1;
    } else if (argc > 3) {
        std::cout << "Too many inputs!" << std::endl;
        return 1;
    }

    const std::string config_path = argv[1];
    const std::string replica_string = argv[2];
    const int replica = stoi(replica_string);

    string targets[] = {"proton", "neutron"};
    string projectiles[] = {"neutrino", "antineutrino"};
    string sf_types[] = {"light", "charm", "top"};

    std::cout << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "Config Path: " << config_path << std::endl;
    std::cout << "Making all splines!" << std::endl;
    std::cout << "REPLICA = " << replica_string << std::endl;
    std::cout << "=============================================" << std::endl << std::endl;

    // Create a new config w/ the filename
    nuxssplmkr::Configuration config = nuxssplmkr::Configuration(config_path);
    config.Populate();  // Populate the SF info
    config.Set_Replica(replica);

    boost::filesystem::path out_folder = "../data/" + config.general.unique_name + "/replica_" + std::to_string(config.pdf.replica);
    if (!boost::filesystem::exists(out_folder)) {  // make the out_folder if it does not exist
        boost::filesystem::create_directories(out_folder);
    }

    std::cout << "Structure function file will be saved to: " << out_folder.string() << std::endl;

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

                sf.InitializeAPFEL();
                sf.BuildSplines(out_folder.string());
            }
        }

    }

    return 0;
}