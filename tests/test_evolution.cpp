#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cmath>
#include "APFEL/APFEL.h"
using namespace std;

int main(int argc, char* argv[]){
    if (argc < 4) {
        std::cout << "Not enough inputs." << std::endl;
        return 0;
    }
    std::string inpdf = argv[1];
    std::string outpdf = argv[2];
    int replica = std::stoi(argv[3]);
    int pto = std::stoi(argv[4]);
    // double xlha[] = {1.1e-9}; //, 5e-2, 1e-1, 2e-1, 3e-1, 4e-1, 5e-1, 6e-1, 7e-1, 8e-1, 9e-1};
    // double Qin = sqrt(2.0);
    
    // APFEL::SetPDFSet("CT18ANNLO");
    APFEL::SetPDFSet(inpdf);
    APFEL::SetReplica(replica);
    APFEL::SetTheory("QCD");
    APFEL::SetQLimits(1.0, 1e4);
    APFEL::SetPerturbativeOrder(pto);
    APFEL::SetPDFEvolution("exactalpha");
    APFEL::SetFastEvolution(false);
    APFEL::SetAlphaEvolution("exact");
    APFEL::SetVFNS();
    APFEL::SetMaxFlavourPDFs(6);
    APFEL::SetMaxFlavourAlpha(6);
    APFEL::SetNumberOfGrids(1);
    APFEL::SetGridParameters(1, 90, 3, 1e-9);
    APFEL::SetPoleMasses(1.3, 4.75, 172.0);
    APFEL::SetAlphaQCDRef(0.118, 91.1870);
    // APFEL::SetLHgridParameters(100, 50, 1e-9, 1e-1, 1, 50, 3, 1e12 );
    double q2min = 1;
    APFEL::SetLHgridParameters(100, 50, 1e-9, 1e-1, 1, 100, q2min, 1e14 );
    // APFEL::SetSmallxResummation(true, "NLL");
    // APFEL::LHAPDFgrid(replica, 1.295, "CT18ANNLO_nf6");
    APFEL::LHAPDFgrid(replica, 1.3, outpdf);
    
    return 0;
}