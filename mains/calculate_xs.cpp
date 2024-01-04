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
    const std::string target = argv[2]; // proton or neutron
    const std::string projectile = argv[3]; // neutrino or antineutrino
    const std::string xs_type = argv[4]; // Which SFs to use total, light, charm, ..

    // Create a new config w/ the filename
    std::cout << config_path << std::endl;
    Configuration config = Configuration(config_path);
    config.Populate();  // Populate the SF info
    
    std::string data_folder = "../data/" + config.sf_info.pdfset + "_" + config.sf_info.mass_scheme + "_pto" + to_string(config.sf_info.perturbative_order);

    config.Set_Target(target);
    config.Set_Projectile(projectile);
    config.Set_SF_Type(xs_type);

    PhysConst* pc = new PhysConst();
    CrossSection* xs = new CrossSection(config);

    // load the three structure function fit files
    string f1 = data_folder + "/F1_" + target + "_" + xs_type + ".fits";
    string f2 = data_folder + "/F2_" + target + "_" + xs_type + ".fits";
    string f3 = data_folder + "/F3_" + target + "_" + xs_type + ".fits";

    // TOTAL XS
    double logemin = 2;
    double logemax = 9;
    int NE = 100;
    double dE = (logemax - logemin) / (NE-1);

    std::ofstream outfile;
    outfile.open(data_folder + "/cross_sections/" + target + "_" + xs_type + ".out");
    
    if (config.sf_info.mass_scheme != "parton") {
        xs->Load_Structure_Functions(f1, f2, f3);
    }

    for (int ei = 0; ei < NE; ei++) {
        double E = pc->GeV * std::pow(10, logemin + ei * dE);
        double _xs = std::log10(xs->TotalXS(E));
        outfile << E << "," << _xs << "\n";
    }
    outfile.close();

    // DOUBLE DIFFERENTIAL XS
    // int Nx = 50;
    // int Ny = 100;

    // double logxmin = -4;
    // double logxmax = 0;
    // double logymin = -4;
    // double logymax = 0;
    // double dx = (logxmax - logxmin) / Nx;
    // double dy = (logymax - logymin) / Ny;
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