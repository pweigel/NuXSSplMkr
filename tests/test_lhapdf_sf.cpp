#include "NuXSSplMkr/configuration.h"
#include "NuXSSplMkr/physconst.h"
#include "NuXSSplMkr/structure_function.h"
#include "NuXSSplMkr/cross_section.h"
#include "NuXSSplMkr/phase_space.h"
#include "NuXSSplMkr/lhapdf_maker.h"
#include "photospline/splinetable.h"
#include "APFEL/APFEL.h"

#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>
#include <chrono>

using namespace nuxssplmkr;

int main(int argc, char* argv[]){
    std::string config_path = argv[1];
    int mode = std::stoi(argv[2]);

    std::string target = "isoscalar";
    std::string projectile = "neutrino";
    std::string sf_type = "light";

    Configuration config = Configuration(config_path);
    config.Populate();

    std::string data_folder = config.general.data_path + "/" + config.general.unique_name + "/replica_" + std::to_string(config.pdf_info.replica);
    config.Set_Target(target);
    config.Set_Projectile(projectile);
    config.Set_SF_Type(sf_type);
    config.Set_Mode(mode);
    LHAPDFMaker lhapdfmaker = LHAPDFMaker(config);
    std::vector<std::string> codes = lhapdfmaker.MakeSet(data_folder);
    lhapdfmaker.MakeInfo(codes);

}