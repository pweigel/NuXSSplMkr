#include "configuration.h"
#include "physconst.h"
#include "structure_function.h"
#include <boost/filesystem.hpp>

int main(int argc, char* argv[]){
    if (argc < 5) {
        std::cout << "Not enough inputs." << std::endl;
        return 0;
    }

    const std::string config_path = argv[1];
    const std::string target = argv[2]; // proton or neutron
    const std::string projectile = argv[3]; // neutrino or antineutrino
    const std::string sf_type = argv[4]; // total, light, charm, ..

    std::cout << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "Config Path: " << config_path << std::endl;
    std::cout << "Target: " << target << std::endl;
    std::cout << "Projectile: " << projectile << std::endl;
    std::cout << "SF Type: " << sf_type << std::endl;
    std::cout << "=============================================" << std::endl << std::endl;

    // Create a new config w/ the filename
    nuxssplmkr::Configuration config = nuxssplmkr::Configuration(config_path);
    config.Populate();  // Populate the SF info

    config.Set_Target(target);
    config.Set_Projectile(projectile);
    config.Set_SF_Type(sf_type);

    nuxssplmkr::StructureFunction sf = nuxssplmkr::StructureFunction(config);
    nuxssplmkr::PhysConst* pc = new nuxssplmkr::PhysConst();

    sf.Set_Lepton_Mass(pc->muon_mass);
    
    std::string _out_folder = "../data/" + config.sf_info.pdfset + "_" + config.sf_info.mass_scheme + "_pto" + to_string(config.sf_info.perturbative_order);
    boost::filesystem::path out_folder = _out_folder;
    if (! boost::filesystem::exists(out_folder)) {
        boost::filesystem::create_directories(out_folder);
    }
    
    sf.InitializeAPFEL();
    sf.BuildSplines(out_folder.string()); // Photospline
    // sf.BuildGrids(out_folder.string()); // Grid file

    return 0;
}