#include "configuration.h"
#include "physconst.h"
#include "structure_function.h"

int main(int argc, char* argv[]){
    const std::string config_path = argv[1];

    // Create a new config w/ the filename
    std::cout << config_path << std::endl;
    nuxssplmkr::Configuration config = nuxssplmkr::Configuration(config_path);
    config.Populate();  // Populate the SF info

    nuxssplmkr::StructureFunction sf = nuxssplmkr::StructureFunction(config);
    std::cout << sf.F2(0., 0.) << std::endl;

    nuxssplmkr::PhysConst* pc = new nuxssplmkr::PhysConst();
    // double enu = 100.;
    // sf.Set_Neutrino_Energy(enu*pc->GeV);
    sf.Set_Lepton_Mass(pc->muon_mass);

    std::string filename_dsdxdy = "test.dat";
    ofstream outfile(filename_dsdxdy.c_str());
    double cm2 = (pc->cm)*(pc->cm);
    for (double logenu=0; logenu<=7.; logenu+=0.05){
        double enu = pow(10, logenu);
        sf.Set_Neutrino_Energy(enu*pc->GeV);
        double sigma = sf.TotalXS();
        outfile << enu << "\t"<< sigma/cm2 << std::endl;

    }
    outfile.close();

    // sf.InitializeAPFEL();

    return 0;
}