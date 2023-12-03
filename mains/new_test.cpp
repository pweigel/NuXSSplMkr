#include "configuration.h"
#include "physconst.h"
#include "structure_function.h"

int main(int argc, char* argv[]){
    const std::string config_path = argv[1];

    // Create a new config w/ the filename
    std::cout << config_path << std::endl;
    nuxssplmkr::Configuration config = nuxssplmkr::Configuration(config_path);
    config.Populate();  // Populate the SF info

    nuxssplmkr::StructureFunction sf = nuxssplmkr::StructureFunction(config);
    std::cout << sf.F2(0., 0.) << std::endl;

    nuxssplmkr::PhysConst* pc = new nuxssplmkr::PhysConst();
    double enu = 100.;
    sf.Set_Neutrino_Energy(enu*pc->GeV);
    sf.Set_Lepton_Mass(enu*pc->muon_mass);

    // sf.InitializeAPFEL();

    return 0;
}