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
    string base_path = "../data/HERAPDF15NLO_EIG_ZM-VFNS_pto1";
    string f1p = base_path + "/F1_proton.fits";
    string f2p = base_path + "/F2_proton.fits";
    string f3p = base_path + "/F3_proton.fits";

    string f1n = base_path + "/F1_neutron.fits";
    string f2n = base_path + "/F2_neutron.fits";
    string f3n = base_path + "/F3_neutron.fits";

    // 'HERAPDF15NLO_EIG_ZM-VFNS_pto1'


    int Nx = 50;
    int Ny = 100;

    double logxmin = -4;
    double logxmax = 0;
    double logymin = -4;
    double logymax = 0;
    double dx = (logxmax - logxmin) / Nx;
    double dy = (logymax - logymin) / Ny;

    double logemin = 2;
    double logemax = 9;
    int NE = 100;
    double dE = (logemax - logemin) / (NE-1);

    std::ofstream proton_outfile;
    proton_outfile.open(base_path + "/cross_sections/proton.out");
    xs->Load_Structure_Functions(f1p, f2p, f3p);
    for (int ei = 0; ei < NE; ei++) {
        double E = pc->GeV * std::pow(10, logemin + ei * dE);
        double _xs = std::log10(xs->TotalXS(E));
        proton_outfile << E << "," << _xs << "\n";
    }
    proton_outfile.close();

    std::ofstream neutron_outfile;
    neutron_outfile.open(base_path + "/cross_sections/neutron.out");
    xs->Load_Structure_Functions(f1n, f2n, f3n);
    for (int ei = 0; ei < NE; ei++) {
        double E = pc->GeV * std::pow(10, logemin + ei * dE);
        double _xs = std::log10(xs->TotalXS(E));
        neutron_outfile << E << "," << _xs << "\n";
    }
    neutron_outfile.close();

    // std::ofstream outfile;
    // outfile.open("test_dsdxdy.out");

    // for (int xi = 0; xi < Nx; xi++) {
    //     double x = logxmin + xi * dx;
    //     for (int yi = 0; yi < Ny; yi++) {
    //         double y = logymin + yi * dy;
    //         // std::cout << std::pow(10, x) << " " << std::pow(10, y) << std::endl;
    //         outfile << xs->ds_dxdy(std::pow(10, x), std::pow(10, y), 100 * pc->GeV);
    //         if (yi != Ny - 1) {
    //             outfile << ",";
    //         }
    //     }
    //     outfile << "\n";
    // }
    // outfile.close();

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