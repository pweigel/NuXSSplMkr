#include "NuXSSplMkr/configuration.h"
#include "NuXSSplMkr/physconst.h"
#include "NuXSSplMkr/structure_function.h"
#include "NuXSSplMkr/spline_maker.h"
#include <boost/filesystem.hpp>

int main(int argc, char* argv[]){
    // if (argc < 5) {
    //     std::cout << "Not enough inputs!" << std::endl;
    //     return 1;
    // } else if (argc > 5) {
    //     std::cout << "Too many inputs!" << std::endl;
    //     return 1;
    // }

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
    std::string data_folder = config.general.data_path + "/" + config.general.unique_name + "/replica_" + std::to_string(config.pdf.replica);
    std::cout << "Loading/saving data to: " << data_folder << std::endl;
    
    // Make the sf folder if it doesn't exist
    boost::filesystem::path out_folder = data_folder;
    if (!boost::filesystem::exists(out_folder)) {
        boost::filesystem::create_directories(out_folder);
    }

    config.Set_Target(target);
    config.Set_Projectile(projectile);
    config.Set_SF_Type(sf_type);

    nuxssplmkr::StructureFunction sf = nuxssplmkr::StructureFunction(config);
    nuxssplmkr::PhysConst* pc = new nuxssplmkr::PhysConst();
    sf.Set_Mode(mode);

    sf.Set_Lepton_Mass(pc->muon_mass);

    if (mode > 1) {
        sf.LoadSplines(data_folder);
    } else {
        sf.InitializeAPFEL();
    }
    sf.BuildGrids(data_folder); // Grid file

    // Spline the structure functions
    nuxssplmkr::SplineMaker splmkr = nuxssplmkr::SplineMaker();
    splmkr.MakeSpline(sf.f1_grid_fn+".grid", sf.f1_grid_fn+".fits");
    splmkr.MakeSpline(sf.f2_grid_fn+".grid", sf.f2_grid_fn+".fits");
    splmkr.MakeSpline(sf.f3_grid_fn+".grid", sf.f3_grid_fn+".fits");

    return 0;
}