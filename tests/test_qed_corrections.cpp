#include "configuration.h"
#include "physconst.h"
#include "structure_function.h"
#include "cross_section.h"
#include "phase_space.h"

#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>

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

    string sf_type = "light";
    string projectile = "neutrino";
    string target = "proton";

    config.Set_SF_Type(sf_type);
    config.Set_Projectile(projectile);
    config.Set_Target(target);

    PhaseSpace ps(config);
    ps.Print();

    std::string data_folder = "../data/" + config.general.unique_name + "/replica_" + std::to_string(config.pdf.replica);

    string f1 = data_folder + "/F1_" + projectile + "_" + target + "_" + sf_type + ".fits";
    string f2 = data_folder + "/F2_" + projectile + "_" + target + "_" + sf_type + ".fits";
    string f3 = data_folder + "/F3_" + projectile + "_" + target + "_" + sf_type + ".fits";

    CrossSection* xs = new CrossSection(config, ps);

    int Ny = 100;
    int Nx = 100;

    double logymin = -7;
    double logymax = 0;
    double dy = (logymax - logymin) / (Ny-1);

    double logxmin = -7;
    double logxmax = 0;
    double dx = (logxmax - logxmin) / (Nx-1);

    std::ofstream outfile;
    outfile.open(data_folder + "/cross_sections/rc_dsdxdy_" + projectile + "_" + target + "_" + sf_type + ".out");
    std::ofstream norc_outfile;
    norc_outfile.open(data_folder + "/cross_sections/no_rc_dsdxdy_" + projectile + "_" + target + "_" + sf_type + ".out");

    xs->Load_Structure_Functions(f1, f2, f3);
    xs->Set_Lepton_Mass(pc->muon_mass);
    std::string base_spline_path = "/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/pweigel/sandbox/src/NuXSSplMkr/scripts/3d_spline_light.fits";
    xs->rc_load_dsdxdy(base_spline_path);

    // Get the energy, inelasticity values and put them in the header
    double E = 100.0 * pc->TeV;
    outfile << "E," << E << std::endl;
    norc_outfile << "E," << E << std::endl;

    std::vector<double> inelasticity_values;
    outfile << "y";
    for (int yi = 0; yi < Ny; yi++) {
        double y = std::pow(10, logymin + yi * dy);
        inelasticity_values.push_back(y);
        outfile << "," << y;
        norc_outfile << "," << y;
    }
    outfile << std::endl;
    norc_outfile << std::endl;

    std::vector<double> x_values;
    outfile << "x";
    for (int xi = 0; xi < Nx; xi++) {
        double x = std::pow(10, logxmin + xi * dx);
        x_values.push_back(x);
        outfile << "," << x;
        norc_outfile << "," << x;
    }
    outfile << std::endl;
    norc_outfile << std::endl;

    for (int yi = 0; yi < Ny; yi++) { // loop over y
        std::cout.flush();
        std::cout << yi << "/" << Ny << "\r";

        double y = inelasticity_values[yi];
        for (int xi = 0; xi < Nx; xi++) { // loop over x
            double x = x_values[xi];
            double _dsdxdy;
            double _rc_dsdxdy;

            bool valid = ps.Validate(E, x, y);
            
            // if the dsdxdy < 1e-38, don't even bother
            if (!valid) {
                _dsdxdy = 0.0;
                _rc_dsdxdy = 0.0;
            } else {
                _dsdxdy = xs->ds_dxdy(E, x, y);
                _rc_dsdxdy = xs->rc_integrate(E, x, y);
                // if (_dsdxdy > 1e-30) {
                    
                // } else{
                //     _rc_dsdxdy = 0.0;
                // }
            }

            outfile << _rc_dsdxdy;
            norc_outfile << _dsdxdy;
            if ( !(xi == Nx - 1) ) {
                outfile << ",";
                norc_outfile << ",";
            }
        }
        outfile << std::endl;
        norc_outfile << std::endl;
    }

    outfile.close();
    norc_outfile.close();

    return 0;
}