#include "configuration.h"
#include "physconst.h"
#include "structure_function.h"
#include <boost/filesystem.hpp>

int main(int argc, char* argv[]){
    if (argc < 5) {
        std::cout << "Not enough inputs!" << std::endl;
        return 1;
    } else if (argc > 5) {
        std::cout << "Too many inputs!" << std::endl;
        return 1;
    }

    const std::string config_path = argv[1];
    const std::string projectile = argv[2]; // neutrino or antineutrino
    const std::string target = argv[3]; // proton or neutron
    const std::string sf_type = argv[4]; // total, light, charm, ..

    std::cout << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "Config Path: " << config_path << std::endl;
    std::cout << "Projectile: " << projectile << std::endl;
    std::cout << "Target: " << target << std::endl;
    std::cout << "SF Type: " << sf_type << std::endl;
    std::cout << "=============================================" << std::endl << std::endl;

    // Create a new config w/ the filename
    nuxssplmkr::Configuration config = nuxssplmkr::Configuration(config_path);
    config.Populate();  // Populate the SF info
    boost::filesystem::path out_folder = "../data/" + config.unique_name;
    if (!boost::filesystem::exists(out_folder)) {  // make the out_folder if it does not exist
        boost::filesystem::create_directories(out_folder);
    }

    config.Set_Target(target);
    config.Set_Projectile(projectile);
    config.Set_SF_Type(sf_type);

    nuxssplmkr::StructureFunction sf = nuxssplmkr::StructureFunction(config);
    nuxssplmkr::PhysConst* pc = new nuxssplmkr::PhysConst();

    sf.Set_Lepton_Mass(pc->muon_mass);

    sf.InitializeAPFEL();
    sf.BuildSplines(out_folder.string()); // Photospline
    // sf.BuildGrids(out_folder.string()); // Grid file

    return 0;
}