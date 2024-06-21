// #include "LHAPDF/LHAPDF.h"
#include "LHAPDF/GridPDF.h"
#include "APFEL/APFEL.h"
#include <deque>
#include <math.h>
#include <iostream>
#include <fstream>

int main(int argc, char* argv[]){

    string pdfset = "CT18ANLO";
    LHAPDF::PDF* pdfs = LHAPDF::mkPDF(pdfset);

    const std::string mass_scheme = argv[1];

    double mc, mb, mt;
    mc = 1.3;
    mb = 4.5;
    mt = 173.0;

    double mc2 = mc*mc;
    double mb2 = mb*mb;
    double mt2 = mt*mt;

    APFEL::SetPDFSet(pdfset);
    APFEL::SetReplica(0);
    APFEL::SetMassScheme(mass_scheme);

    APFEL::SetPerturbativeOrder(2);

    APFEL::SetPoleMasses(mc, mb, mt);
    APFEL::SetQLimits(1.0, 1e6);
    APFEL::EnableDampingFONLL(true);


    APFEL::SetNumberOfGrids(3);
    APFEL::SetGridParameters(1, 60, 3, 1e-11);
    APFEL::SetGridParameters(2, 40, 5, 1e-1);
    APFEL::SetGridParameters(3, 30, 5, 8e-1);
    LHAPDF::PDF* pdf = LHAPDF::mkPDF(pdfset, 0);

    APFEL::SetAlphaQCDRef(pdf->alphasQ(91.1876), 91.1876);

    double Vud = 0.9743;
    double Vus = 0.2254;
    double Vub = 0.0036;
    double Vcd = 0.2252;
    double Vcs = 0.9734;
    double Vcb = 0.0414;
    double Vtd = 0.0089;
    double Vts = 0.0405;
    double Vtb = 0.9991;

    Vud = 1.;
    Vus = 0.;
    Vub = 0.;
    Vcd = 0.;
    Vcs = 1.;
    Vcb = 0.;
    Vtd = 0.;
    Vts = 0.;
    Vtb = 1.;

    APFEL::SetMaxFlavourPDFs(6);
    APFEL::SetMaxFlavourAlpha(6);
    APFEL::SetCKM(Vud, Vus, Vub,
                  Vcd, Vcs, Vcb,
                  Vtd, Vts, Vtb);
    APFEL::SetProcessDIS("CC");

    // APFEL::SetProjectileDIS("antineutrino");
    APFEL::SetProjectileDIS("neutrino");
    APFEL::SetTargetDIS("proton");

    // Initializes integrals on the grids
    APFEL::InitializeAPFEL_DIS();
    double x = 0.1;
    APFEL::ComputeStructureFunctionsAPFEL(4, 4);


    std::cout << APFEL::F2total(x) << ", " << APFEL::FLtotal(x) << ", " << APFEL::F2total(x) - APFEL::FLtotal(x) << std::endl;

}