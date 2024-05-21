#include "NuXSSplMkr/configuration.h"
#include "NuXSSplMkr/physconst.h"
#include "NuXSSplMkr/structure_function.h"
#include "NuXSSplMkr/spline_maker.h"
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

    int replica = 0;

    // Create a new config w/ the filename
    std::cout << config_path << std::endl;
    nuxssplmkr::Configuration config = nuxssplmkr::Configuration(config_path);
    config.Populate();
    config.Set_Replica(replica);
    std::string data_folder = "../data/" + config.general.unique_name + "/replica_" + std::to_string(config.pdf.replica);
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

    sf.Set_Lepton_Mass(pc->muon_mass);

    sf.InitializeAPFEL();
    // sf.BuildSplines(out_folder.string()); // Photospline
    sf.BuildGrids(out_folder.string()); // Grid file

    // nuxssplmkr::SplineMaker splmkr = nuxssplmkr::SplineMaker(config);
    // string gridinfile = data_folder + "/F1_neutrino_proton_light.grid";
    // string fitsoutfile = data_folder + "/F1_neutrino_proton_light.fits";
    // splmkr.MakeSpline(gridinfile, fitsoutfile);

    return 0;
}