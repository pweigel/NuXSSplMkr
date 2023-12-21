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
    nuxssplmkr::PhysConst* pc = new nuxssplmkr::PhysConst();

    sf.Set_Lepton_Mass(pc->muon_mass);

    std::string sf_grids_fname = "test_sf_grid.dat";
    // x Q2 F1 F2 xF3
    ofstream outfile(sf_grids_fname.c_str());
    double cm2 = (pc->cm)*(pc->cm);

    // for (double logx=-7; logx<0.; logx+=0.1) {
    //     for (double logQ2=0; logQ2<3; logQ2+=0.1) {
    //         double x = pow(10, logx);
    //         double Q2 = pow(10, logQ2);
    //         double F1 = sf.F1(x, Q2);
    //         double F2 = sf.F2(x, Q2);
    //         double xF3 = sf.xF3(x, Q2);
    //         outfile << x << " " << Q2 << " " << F1 << " " << F2 << " " << xF3 << std::endl;
    //     }
    // }

    // sf.Set_Use_APFEL_LO(false);
    sf.InitializeAPFEL();
    for (double logx=-4; logx<0.; logx+=0.01) {
        double x = pow(10, logx);
        double Q2 = 4.;
        
        sf.Set_Q_APFEL(std::sqrt(Q2));

        double F1 = sf.F1(x, Q2);
        double F2 = sf.F2(x, Q2);
        double xF3 = sf.xF3(x, Q2);
        outfile << x << " " << Q2 << " " << F1 << " " << F2 << " " << xF3 << std::endl;
    }

    // for (double logenu=0; logenu<=7.; logenu+=0.05){
    //     double enu = pow(10, logenu);
    //     sf.Set_Neutrino_Energy(enu*pc->GeV);
    //     double sigma = sf.TotalXS();
    //     outfile << enu << "\t"<< sigma/cm2 << std::endl;

    // }
    outfile.close();

    std::cout << "Done!" << std::endl;



    return 0;
}