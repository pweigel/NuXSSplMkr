#include "configuration.h"
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

    nuxssplmkr::Configuration config = nuxssplmkr::Configuration(config_path);
    
    config.Populate();  // Populate the SF info

    config.Set_Target(target);
    config.Set_Projectile(projectile);
    config.Set_SF_Type(sf_type);

    nuxssplmkr::StructureFunction sf = nuxssplmkr::StructureFunction(config);
    nuxssplmkr::PhysConst* pc = new nuxssplmkr::PhysConst();

    sf.Set_Lepton_Mass(pc->muon_mass);

    double Q2 = 2.0;

    int Nx = 100;
    double logxmin = -3;
    double logxmax = 0;
    double dx = (logxmax - logxmin) / (Nx-1);
    sf.InitializeAPFEL();
    sf.Set_Q_APFEL(std::sqrt(Q2));
    std::ofstream outfile;
    outfile.open("../test_data/TMC_"+target+".out");
    for (int i = 0; i < Nx; i++) {
        double x = pow(10, logxmin + i*dx);
        std::cout << x << ",";
        double f2 = sf.F2(x, Q2);
        std::cout << f2 << ",";
        double xf3 = sf.xF3(x, Q2);
        std::cout << xf3 << ",";
        double f2_tmc = sf.F2_TMC(x, Q2);
        std::cout << f2_tmc << ",";
        double xf3_tmc = sf.xF3_TMC(x, Q2);
        std::cout << xf3_tmc << std::endl;
        outfile << x << "," << f2 << "," << xf3 << "," << f2_tmc << "," << xf3_tmc << "\n";
    }

    // sf.BuildSplines(out_folder.string()); // Photospline
    // sf.BuildGrids(out_folder.string()); // Grid file

    return 0;
}