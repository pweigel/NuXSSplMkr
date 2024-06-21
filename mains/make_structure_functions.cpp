#include "NuXSSplMkr/configuration.h"
#include "NuXSSplMkr/physconst.h"
#include "NuXSSplMkr/structure_function.h"
// #include "NuXSSplMkr/spline_maker.h"
#include "NuXSSplMkr/lhapdf_maker.h"

#include <boost/filesystem.hpp>

int main(int argc, char* argv[]){
    if (argc != 7) {
        std::cout << "Not enough/too many inputs!" << std::endl;
        std::cout << "Usage: make_structure_functions CONFIG PROJECTILE TARGET TYPE MODE REPLICA" << std::endl;
        return 1;
    }


    const std::string config_path = argv[1];
    const std::string projectile = argv[2]; // neutrino or antineutrino
    const std::string target = argv[3]; // proton or neutron
    const std::string sf_type = argv[4]; // total, light, charm, ..
    const unsigned int mode = std::stoi(argv[5]);
    const unsigned int replica = std::stoi(argv[6]);

    std::cout << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "Config Path: " << config_path << std::endl;
    std::cout << "Projectile: " << projectile << std::endl;
    std::cout << "Target: " << target << std::endl;
    std::cout << "SF Type: " << sf_type << std::endl;
    std::cout << "Mode: " << mode << std::endl;
    std::cout << "Replica: " << replica << std::endl;
    std::cout << "=============================================" << std::endl << std::endl;

    // Create a new config w/ the filename
    std::cout << config_path << std::endl;
    nuxssplmkr::Configuration config = nuxssplmkr::Configuration(config_path);
    config.Populate();
    config.Set_Replica(replica);
    std::string data_folder = config.general.data_path + "/" + config.general.unique_name + "/replica_" + std::to_string(config.pdf_info.replica);
    std::cout << "Loading/saving data to: " << data_folder << std::endl;
    
    // Make the sf folder if it doesn't exist
    boost::filesystem::path out_folder = data_folder;
    if (!boost::filesystem::exists(out_folder)) {
        boost::filesystem::create_directories(out_folder);
    }

    nuxssplmkr::PhysConst* pc = new nuxssplmkr::PhysConst();

    config.Set_Target(target);
    config.Set_Projectile(projectile);
    config.Set_SF_Type(sf_type);
    config.Set_Lepton_Mass(pc->muon_mass);
    config.Set_Mode(mode);
    
    nuxssplmkr::StructureFunction sf = nuxssplmkr::StructureFunction(config);

    if (mode == 1) {
        sf.InitializeAPFEL();
    }

    sf.BuildGrids(data_folder);

    return 0;
}