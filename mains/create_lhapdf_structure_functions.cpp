#include "NuXSSplMkr/configuration.h"
#include "NuXSSplMkr/lhapdf_maker.h"

using namespace nuxssplmkr;

int main(int argc, char* argv[]){
    if (argc != 4) {
        std::cout << "Not enough/too many inputs!" << std::endl;
        std::cout << "Usage: create_lhapdf_structure_functions CONFIG MODE REPLICA" << std::endl;
        return 1;
    }

    std::string config_path = argv[1];
    const int mode = std::stoi(argv[2]);
    const int replica = std::stoi(argv[3]);

    Configuration config = Configuration(config_path);
    config.Populate();
    config.Set_Replica(replica);

    std::string data_folder = config.general.data_path + "/" + config.general.unique_name + "/replica_" + std::to_string(config.pdf_info.replica);
    config.Set_Mode(mode);
    LHAPDFMaker lhapdfmaker = LHAPDFMaker(config);
    std::vector<std::string> codes = lhapdfmaker.MakeSet(data_folder);

    // Make .info file if it doesn't exist (lhapdfmaker will check if it does)
    lhapdfmaker.MakeInfo(codes);
}