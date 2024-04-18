#include "configuration.h"
#include "physconst.h"
#include "structure_function.h"
#include "cross_section.h"
#include "phase_space.h"
#include "photospline/splinetable.h"

#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>
#include <chrono>

using namespace nuxssplmkr;

int main(int argc, char* argv[]){
    const std::string config_path = argv[1]; // Path to .json file containing configuration info

    std::cout << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "Config Path: " << config_path << std::endl;
    std::cout << "=============================================" << std::endl << std::endl;

    PhysConst* pc = new PhysConst();

    Configuration config = Configuration(config_path);
    config.Populate();
    config.Set_Replica(config.pdf.replica);

    string sf_type = "total";
    string projectile = "neutrino";
    string target = "proton";

    config.Set_SF_Type(sf_type);
    config.Set_Projectile(projectile);
    config.Set_Target(target);

    PhaseSpace ps(config);
    ps.Print();

    std::string data_folder = "../data/" + config.general.unique_name + "/replica_" + std::to_string(config.pdf.replica);

    // string f1 = data_folder + "/F1_" + projectile + "_" + target + "_" + sf_type + ".fits";
    // string f2 = data_folder + "/F2_" + projectile + "_" + target + "_" + sf_type + ".fits";
    // string f3 = data_folder + "/F3_" + projectile + "_" + target + "_" + sf_type + ".fits";

    CrossSection* xs = new CrossSection(config, ps);

    int NE = 80;
    int Ny = 90;
    int Nx = 90;

    double logemin = 2;
    double logemax = 9;
    double dE = (logemax - logemin) / (NE-1);

    double logymin = -9;
    double logymax = -1;
    double dy = (logymax - logymin) / ( (Ny-20) - 1);

    double logxmin = -9;
    double logxmax = 0;
    double dx = (logxmax - logxmin) / (Nx-1);

    std::ofstream outfile;
    // outfile.open(data_folder + "/cross_sections/rc_dsdxdy_" + projectile + "_" + target + "_" + sf_type + ".out");
    outfile.open("dsdxdy_nu_CC_iso_corrected_logspaced.out");
    // std::ofstream norc_outfile;
    // norc_outfile.open(data_folder + "/cross_sections/no_rc_dsdxdy_" + projectile + "_" + target + "_" + sf_type + ".out");

    // xs->Load_Structure_Functions(f1, f2, f3);
    xs->Set_Lepton_Mass(pc->muon_mass);
    // std::string base_spline_path = "/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/pweigel/sandbox/src/NuXSSplMkr/scripts/3d_spline_light.fits";
    std::string base_spline_path = "/n/home06/pweigel/utils/xs_iso/dsdxdy_nu_CC_iso.fits";
    xs->rc_load_dsdxdy(base_spline_path);

    // Get the energy, inelasticity values and put them in the header
    std::vector<double> energy_values;
    outfile << "E";
    for (int ei = 0; ei < NE; ei++) {
        double E = pc->GeV * std::pow(10, logemin + ei * dE);
        energy_values.push_back(E);
        outfile << "," << E;
    }
    outfile << std::endl;

    /*
    Using Ny = 90, use 70 for the log spacing and 20 for linear
    
    */
    std::vector<double> inelasticity_values;
    outfile << "y";
    for (int yi = 0; yi < Ny - 20; yi++) {
        double y = std::pow(10, logymin + yi * dy);
        inelasticity_values.push_back(y);
        outfile << "," << y;
        // norc_outfile << "," << y;
        // std::cout << y << std::endl;
    }

    for (int yi = 1; yi < 20; yi++) {
        double y = std::pow(10, logymin + dy * (Ny-20-1)) + yi * (1.0 - 0.1) / (20);
        inelasticity_values.push_back(y);
        outfile << "," << y;
        // std::cout << y << std::endl;
    }
    inelasticity_values.push_back(0.98);
    outfile << "," << "0.98";
    std::cout << "size of y-values = " << inelasticity_values.size() << std::endl;

    outfile << std::endl;
    // norc_outfile << std::endl;

    std::vector<double> x_values;
    outfile << "x";
    for (int xi = 0; xi < Nx-1; xi++) {
        double x = std::pow(10, logxmin + xi * dx);
        x_values.push_back(x);
        outfile << "," << x;
        // norc_outfile << "," << x;
    }
    x_values.push_back(0.99);
    outfile << "," << "0.99";
    outfile << std::endl;
    // norc_outfile << std::endl;

    auto t0 = std::chrono::high_resolution_clock::now();
    double delta_time = 0.;
    for (int ei = 0; ei < NE; ei++) { // loop over E

        double E = energy_values[ei];
        std::cout.flush();
        std::cout << ei << "/" << NE << ", dt=" << delta_time << " seconds." << "\r";

        auto t1 = std::chrono::high_resolution_clock::now();

        for (int yi = 0; yi < Ny; yi++) { // loop over y
            double y = inelasticity_values[yi];
            for (int xi = 0; xi < Nx; xi++) { // loop over x
                double x = x_values[xi];
                // double _dsdxdy;
                double _rc_dsdxdy;
                double xs_val;

                bool valid = ps.Validate(E, x, y);
                
                if (!valid) {
                    // _dsdxdy = 0.0;
                    _rc_dsdxdy = 0.0;
                    xs_val = 0.0;
                } else {
                    // _dsdxdy = xs->ds_dxdy(E, x, y);
                    _rc_dsdxdy = xs->rc_integrate(E, x, y);
                    std::array<double, 3> pt{{std::log10(E / pc->GeV), std::log10(x), std::log10(y)}};
                    std::array<int, 3> xs_splc;
                    xs->rc_dsdxdy.searchcenters(pt.data(), xs_splc.data());
                    xs_val = pow(10.0, xs->rc_dsdxdy.ndsplineeval(pt.data(), xs_splc.data(), 0));
                }

                outfile << xs_val + _rc_dsdxdy;
                // norc_outfile << _dsdxdy;
                if ( !(xi == Nx - 1) ) {
                    outfile << ",";
                    // norc_outfile << ",";
                }
            }
            outfile << std::endl;
            // norc_outfile << std::endl;
            auto t2 = std::chrono::high_resolution_clock::now();

            delta_time = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() / 1e3;
        }
    }

    outfile.close();
    // norc_outfile.close();
    auto t3 = std::chrono::high_resolution_clock::now();
    delta_time = std::chrono::duration_cast<std::chrono::milliseconds>(t3-t0).count() / 1e3;
    std::cout << "TOTAL TIME = " << delta_time << " seconds." << std::endl;
    return 0;
}