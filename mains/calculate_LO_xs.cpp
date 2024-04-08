#include "configuration.h"
#include "physconst.h"
#include "structure_function.h"
#include <boost/filesystem.hpp>

int main(int argc, char* argv[]){
    if (argc < 2) {
        std::cout << "Not enough inputs." << std::endl;
        return 0;
    }
    const std::string config_path = argv[1];

    // Create a new config w/ the filename
    std::cout << config_path << std::endl;
    nuxssplmkr::Configuration config = nuxssplmkr::Configuration(config_path);
    config.Populate();  // Populate the SF info

    nuxssplmkr::StructureFunction sf = nuxssplmkr::StructureFunction(config);
    nuxssplmkr::PhysConst* pc = new nuxssplmkr::PhysConst();

    sf.Set_Lepton_Mass(pc->muon_mass);
    
    std::string _out_folder = "../data/" + config.sf_info.pdfset + "_" + config.sf_info.mass_scheme + "_pto" + to_string(config.sf_info.perturbative_order);
    boost::filesystem::path out_folder = _out_folder;
    if (! boost::filesystem::exists(out_folder)) {
        boost::filesystem::create_directories(out_folder);
    }

    int Nx = 100;
    int Ny = 100;
    double logxmin = -4;
    double logxmax = 0;
    double logymin = -4;
    double logymax = 0;
    double dx = (logxmax - logxmin) / Nx;
    double dy = (logymax - logymin) / Ny;

    std::ofstream outfile;
    outfile.open("LO_dsdxdy.out");
    // outfile.open("test_totalxs.out");
    double cm2 = (pc->cm) * (pc->cm);
    for (int xi = 0; xi < Nx; xi++) {
        double x = logxmin + xi * dx;
        for (int yi = 0; yi < Ny; yi++) {
            double y = logymin + yi * dy;
            // std::cout << std::pow(10, x) << " " << std::pow(10, y) << std::endl;
            outfile << sf.ds_dxdy_LO(std::pow(10, x), std::pow(10, y), 100 * pc->GeV);

            if (yi != Ny - 1) {
                outfile << ",";
            }
        }
        outfile << "\n";
    }
    outfile.close();

    return 0;
}