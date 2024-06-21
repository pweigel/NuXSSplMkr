#include "NuXSSplMkr/configuration.h"
#include "NuXSSplMkr/physconst.h"
#include "NuXSSplMkr/structure_function.h"
#include "NuXSSplMkr/cross_section.h"
#include "NuXSSplMkr/phase_space.h"
#include "photospline/splinetable.h"
#include "APFEL/APFEL.h"

#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>
#include <chrono>

using namespace nuxssplmkr;

int main(int argc, char* argv[]){
    const std::string config_path = argv[1]; // Path to .json file containing configuration info
    const std::string projectile = argv[2]; // neutrino or antineutrino
    const unsigned int replica = 0;

    PhysConst* pc = new PhysConst();
    Configuration config = Configuration(config_path);

    config.Populate();
    config.Set_Replica(replica);
    std::string data_folder = config.general.data_path + "/" + config.general.unique_name + "/replica_" + std::to_string(replica);

    config.Set_Projectile("neutrino");
    config.Set_SF_Type("light");
    config.Set_Target("proton");
    PhaseSpace ps = PhaseSpace(config);

    CrossSection* xs = new CrossSection(config, ps);
    xs->Load_Structure_Functions(data_folder);
    xs->Set_Lepton_Mass(pc->muon_mass);

    double energies[] = {1e3*pc->GeV, 1e4*pc->GeV, 1e5*pc->GeV, 1e6*pc->GeV};

    for (int i = 0; i < 4; i++) {
        double E = energies[i];
    
        // double E = 1e4 * pc->GeV;

        double xs_a = xs->TotalXS(E) /1e-36;
        double xs_b = xs->TotalXS_xQ2(E) /1e-36;
        double xs_c = xs->AlternativeTotalXS(E) /1e-36;
        // std::cout << E << ": " << xs_b/xs_a << ", " << xs_c/xs_a << std::endl;
        std::cout << E << ": " << xs_a << ", " << xs_b << ", " << xs_c << std::endl;
    }
    return 0;
}