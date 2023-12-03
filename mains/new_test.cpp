#include "configuration.h"
#include "structure_function.h"

int main(int argc, char* argv[]){
    const std::string config_path = argv[1];

    // Create a new config w/ the filename
    std::cout << config_path << std::endl;
    nuxssplmkr::Configuration config = nuxssplmkr::Configuration(config_path);
    config.Populate();  // Populate the SF info

    nuxssplmkr::StructureFunction sf = nuxssplmkr::StructureFunction(config);
    std::cout << sf.F2(0., 0.) << std::endl;

    return 0;
}