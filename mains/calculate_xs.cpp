#include "configuration.h"
#include "physconst.h"
#include "structure_function.h"
#include "cross_section.h"
#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>

using namespace nuxssplmkr;

int main(int argc, char* argv[]){
    if (argc < 2) {
        std::cout << "Not enough inputs." << std::endl;
        return 0;
    }
    const std::string config_path = argv[1];

    // Create a new config w/ the filename
    std::cout << config_path << std::endl;
    Configuration config = Configuration(config_path);
    config.Populate();  // Populate the SF info

    PhysConst* pc = new PhysConst();
    CrossSection* xs = new CrossSection(config);

    // load the three structure function fit files
    // string f1 = "../data/CT18NLO_NF6_FONLL-B_pto1/F1.fits";
    // string f2 = "../data/CT18NLO_NF6_FONLL-B_pto1/F2.fits";
    // string f3 = "../data/CT18NLO_NF6_FONLL-B_pto1/F3.fits";
    // string f1 = "../data/CT18NLO_NF6_ZM-VFNS_pto1/F1.fits";
    // string f2 = "../data/CT18NLO_NF6_ZM-VFNS_pto1/F2.fits";
    // string f3 = "../data/CT18NLO_NF6_ZM-VFNS_pto1/F3.fits";
    string f1 = "../data/HERAPDF15NLO_EIG_ZM-VFNS_pto1/F1.fits";
    string f2 = "../data/HERAPDF15NLO_EIG_ZM-VFNS_pto1/F2.fits";
    string f3 = "../data/HERAPDF15NLO_EIG_ZM-VFNS_pto1/F3.fits";
    // 'HERAPDF15NLO_EIG_ZM-VFNS_pto1'
    xs->Load_Structure_Functions(f1, f2, f3);


    std::ofstream outfile;
    outfile.open("test_dsdxdy.out");
    int Nx = 100;
    int Ny = 100;
    double logxmin = -4;
    double logxmax = 0;
    double logymin = -4;
    double logymax = 0;
    double dx = (logxmax - logxmin) / Nx;
    double dy = (logymax - logymin) / Ny;
    double E = 1000.0 * pc->GeV;
    std::cout << xs->ds_dxdy(0.1, 0.1, 100.0 * pc->GeV) << std::endl;

    std::cout << "Calculating total cross section!" << std::endl;
    double total_xs = xs->TotalXS(E);
    std::cout << total_xs << std::endl;

    for (int xi = 0; xi < Nx; xi++) {
        double x = logxmin + xi * dx;
        for (int yi = 0; yi < Ny; yi++) {
            double y = logymin + yi * dy;
            // std::cout << std::pow(10, x) << " " << std::pow(10, y) << std::endl;
            outfile << xs->ds_dxdy(std::pow(10, x), std::pow(10, y), E);
            if (yi != Ny - 1) {
                outfile << ",";
            }
        }
        outfile << "\n";
    }
    outfile.close();

    // std::cout << xs->ds_dxdy(0.1, 0.1, 100.0 * pc->GeV) << std::endl;
    // std::cout << xs->ds_dxdy(0.1, 0.1, 1000.0 * pc->GeV) << std::endl;
    // std::cout << xs->ds_dxdy(0.1, 0.1, 10000.0 * pc->GeV) << std::endl;
    // std::string _out_folder = "../data/" + config.sf_info.pdfset + "_" + config.sf_info.mass_scheme + "_pto" + to_string(config.sf_info.perturbative_order);
    // boost::filesystem::path out_folder = _out_folder;
    // if (! boost::filesystem::exists(out_folder)) {
    //     boost::filesystem::create_directories(out_folder);
    // }
    
    // sf.InitializeAPFEL();
    // sf.BuildGrids(out_folder.string());

    return 0;
}