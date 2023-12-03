#include "structure_function.h"

namespace nuxssplmkr {

StructureFunction::StructureFunction(Configuration &_config)
  : config(_config)
{
    pc = new PhysConst();
    sf_info = config.sf_info;

    // TODO: Get from config file
    GF2 = SQ(pc->GF);
    M_iso = 0.5*(pc->proton_mass + pc->neutron_mass);
    Mw2 = SQ(pc->Wboson_mass);
    Mz2 = SQ(pc->Zboson_mass);

    // Calculate fundamental constants
    s_w = sf_info.Sin2ThW;
    Lu2 = ( 1. - (4./3.)*s_w) * ( 1. - (4./3.)*s_w);
    Ld2 = (-1. + (2./3.)*s_w) * (-1. + (2./3.)*s_w);
    Ru2 = (    - (4./3.)*s_w) * (    - (4./3.)*s_w);
    Rd2 = (      (2./3.)*s_w) * (      (2./3.)*s_w);

}

void StructureFunction::InitializeAPFEL() {
    //APFEL::SetMassScheme("ZM-VFNS");
    APFEL::SetPDFSet(sf_info.pdfset);
    APFEL::SetReplica(sf_info.replica);
    APFEL::SetMassScheme(sf_info.mass_scheme);

    APFEL::SetProcessDIS(sf_info.DIS_process);
    APFEL::SetQLimits(std::sqrt(sf_info.Q2min), std::sqrt(sf_info.Q2max));
    //APFEL::SetPolarizationDIS(0);

    APFEL::SetProjectileDIS(sf_info.projectile);
    APFEL::SetTargetDIS(sf_info.target);
    //APFEL::EnableTargetMassCorrections(false);
    //APFEL::EnableDampingFONLL(true);
    //APFEL::SetFastEvolution(true);
    //APFEL::LockGrids(true);
    //APFEL::EnableEvolutionOperator(true);
    //APFEL::SetFFNS(3);
    //APFEL::SetTheory("QUniD");
    //APFEL::SetTheory("QED");
    //APFEL::SetTheory("QUniD");
    //APFEL::EnableLeptonEvolution(true);
    //APFEL::SetTauMass(1e10);
    //APFEL::SetPDFEvolution("exactalpha");
    APFEL::SetNumberOfGrids(3);
    APFEL::SetGridParameters(1, 30, 3, sf_info.xmin);
    APFEL::SetGridParameters(2, 30, 5, 2e-1);
    APFEL::SetGridParameters(3, 30, 5, 8e-1);
    APFEL::SetPerturbativeOrder(sf_info.perturbative_order);
    APFEL::SetAlphaQCDRef(sf_info.pdf->alphasQ(sf_info.MassZ), sf_info.MassZ);
    //APFEL::SetAlphaEvolution("expanded");
    //APFEL::SetPDFEvolution("expandalpha");
    APFEL::SetPoleMasses(sf_info.pdf_quark_masses[4], sf_info.pdf_quark_masses[5], sf_info.pdf_quark_masses[6]);
    APFEL::SetMaxFlavourPDFs(6);
    APFEL::SetMaxFlavourAlpha(6);
    APFEL::SetCKM(sf_info.Vud, sf_info.Vus, sf_info.Vub,
                  sf_info.Vcd, sf_info.Vcs, sf_info.Vcb,
                  sf_info.Vtd, sf_info.Vts, sf_info.Vtb);

    // Initializes integrals on the grids
    APFEL::InitializeAPFEL_DIS();
}

double StructureFunction::F1(double x, double Q2) {
    // LO for now
    return F2(x, Q2) / (2. * x);
}

double StructureFunction::F2(double x, double Q2) {
    // LO for now
    auto s = PDFExtract(x, Q2);
    return F2_LO(s);
}

double StructureFunction::F2_LO(map<int, double>& xq_arr) {
    double k = 0.;

    map<int,double> F2coef;
    F2coef[1]  = 1.;
    F2coef[-1] = 1.;
    F2coef[2]  = 1.;
    F2coef[-2] = 1.;
    F2coef[3]  = 2.;
    F2coef[-3] = 0.;
    F2coef[4]  = 0.;
    F2coef[-4] = 2.;
    F2coef[5]  = 2.;
    F2coef[-5] = 0.;
    F2coef[21] = 0.;

    // mean value
    for( int p : partons ) {
        k += F2coef[p] * xq_arr[p];
    }	
    
    return k;
}

double StructureFunction::xF3(double x, double Q2) {
    // only true at LO
    auto s = PDFExtract(x, Q2);
    return xF3_LO(s);
}

double StructureFunction::xF3_LO(map<int, double>& xq_arr) {
    double k=0.;

    // xF3 coeficients
    map<int,double> F3coef;
    F3coef[1]  = 1.;
    F3coef[-1] = -1.;
    F3coef[2]  = 1.;
    F3coef[-2] = -1.;
    F3coef[3]  = 2.;
    F3coef[-3] = 0.;
    F3coef[4]  = 0.;
    F3coef[-4] = -2.;
    F3coef[5]  = 2.;
    F3coef[-5] = 0.;
    F3coef[21] = 0.;

    // mean value
    for( int p : partons ) {
        k += F3coef[p]*xq_arr[p];
    }	

    return k;
}

double StructureFunction::F3(double x, double Q2) {
    return xF3(x, Q2) / x;
}

std::map<int,double> StructureFunction::PDFExtract(double x, double Q2){
    LHAPDF::GridPDF* grid_central = dynamic_cast<LHAPDF::GridPDF*>(sf_info.pdf);
    string xt = "nearest";
    grid_central -> setExtrapolator(xt);

    std::map<int,double> xq_arr;
    for ( int p : partons ){
      xq_arr[p] = grid_central -> xfxQ2(p, x, Q2/SQ(pc->GeV));
    }
    return xq_arr;
}

double StructureFunction::CrossSection(double x, double Q2) {
    // To be replaced with a different object

    return 0.;
}

double StructureFunction::Evaluate(double Q2, double x, double y){
    // only evaluates central values
    double q = sqrt(Q2)/pc->GeV;

    LHAPDF::GridPDF* grid_central = dynamic_cast<LHAPDF::GridPDF*>(sf_info.pdf);
    string xt = "nearest";
    grid_central -> setExtrapolator(xt);

    map<int,double> xq_arr;
    for ( int p : partons ){
      xq_arr[p] = grid_central -> xfxQ(p,x,q);
    }
    
    return SigR_Nu_LO(x, y, xq_arr);
}

double StructureFunction::SigR_Nu_LO(double x, double y, map<int,double> xq_arr){
	double k = 0.;
    d_lepton = SQ(M_lepton)/(2.*M_iso*ENU);
    double CP_factor = 1.;
	double y_p = (1. - d_lepton / x) + (1.- d_lepton/x - y) * (1. - y);
	double y_m = (1. - d_lepton / x) - (1.- d_lepton/x - y) * (1. - y);
	double a = y_p + CP_factor*y_m;
	double b = y_p - CP_factor*y_m;

    map<int,double> SigRcoef;

    // Coefficients for CC
	SigRcoef[1]  =    a ;
	SigRcoef[-1] =    b ;
	SigRcoef[2]  =    a ;
	SigRcoef[-2] =    b ;
	SigRcoef[3]  = 2.*a ;
	SigRcoef[-3] =   0. ;
	SigRcoef[4]  =   0. ;
	SigRcoef[-4] = 2.*b ;
	SigRcoef[5]  = 2.*a ;
	SigRcoef[-5] =   0. ;
	SigRcoef[21] =   0. ;

    if (CP_factor < 0 ){
        //fixes for antineutrinos
        SigRcoef[3]   = 0. ;
        SigRcoef[-3]  = 2.*b ;
        SigRcoef[4]   = 2.*a ;
        SigRcoef[-4]  = 0. ;
        SigRcoef[5]   = 0. ;
        SigRcoef[-5]  = 2.*b ;
    }

    for( int p : partons ) {
        k += SigRcoef[p]*xq_arr[p];
    }

    return k;
}

void StructureFunction::Set_Lepton_Mass(double m) {
    M_lepton = m;
}

void StructureFunction::Set_Neutrino_Energy(double E) {
    ENU = E;
}

}