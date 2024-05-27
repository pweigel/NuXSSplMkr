#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cmath>
#include "APFEL/APFEL.h"
using namespace std;

int main()
{
    double xlha[] = {1.1e-9}; //, 5e-2, 1e-1, 2e-1, 3e-1, 4e-1, 5e-1, 6e-1, 7e-1, 8e-1, 9e-1};
    double Qin = sqrt(2.0);
    
    APFEL::SetPDFSet("CT18ANNLO");
    APFEL::SetReplica(0);
    APFEL::SetTheory("QCD");
    APFEL::SetPerturbativeOrder(2);
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
    APFEL::SetLHgridParameters(100, 50, 1e-9, 1e-1, 1, 50, 1.295, 1e10 );
    APFEL::LHAPDFgrid(0, 1.295, "aasdf");
    
    return 0;
    APFEL::InitializeAPFEL();

    double Q0, Q, eps = 1e-10;
    cout << "Enter initial and final scale in GeV" << endl;
    cin >> Q0 >> Q;

    // Load evolution
    Q0 = Q0 - eps;
    APFEL::EvolveAPFEL(Q0, Q);

    std::cout << std::scientific << std::setprecision(5) << std::endl;
    cout << "   x   " 
        << setw(11) << "     tbar    " 
        << setw(11) << "     bbar    " 
        << setw(11) << "     cbar    " 
        << setw(11) << "     sbar    " 
        << setw(11) << "     ubar    " 
        << setw(11) << "     dbar    " 
        << setw(11) << "       d     " 
        << setw(11) << "       u     " 
        << setw(11) << "       s     " 
        << setw(11) << "       c     " 
        << setw(11) << "       b     " 
        << setw(11) << "       t     " 
        << setw(11) << "       g     " 
        << endl;

    cout << scientific;
    for (int i = 1; i < 2; i++)
        cout << setprecision(1) << xlha[i] << "\t" << setprecision(4) 
        << setw(11) << APFEL::xPDF(-6,xlha[i]) << "  "
        << setw(11) << APFEL::xPDF(-5,xlha[i]) << "  "
        << setw(11) << APFEL::xPDF(-4,xlha[i]) << "  "
        << setw(11) << APFEL::xPDF(-3,xlha[i]) << "  "
        << setw(11) << APFEL::xPDF(-2,xlha[i]) << "  "
        << setw(11) << APFEL::xPDF(-1,xlha[i]) << "  "
        << setw(11) << APFEL::xPDF(1,xlha[i]) << "  "
        << setw(11) << APFEL::xPDF(2,xlha[i]) << "  "
        << setw(11) << APFEL::xPDF(3,xlha[i]) << "  "
        << setw(11) << APFEL::xPDF(4,xlha[i]) << "  "
        << setw(11) << APFEL::xPDF(5,xlha[i]) << "  "
        << setw(11) << APFEL::xPDF(6,xlha[i]) << "  "
        << setw(11) << APFEL::xPDF(0,xlha[i]) << "  "
        << endl;
    cout << "      " << endl;
    return 0;
}