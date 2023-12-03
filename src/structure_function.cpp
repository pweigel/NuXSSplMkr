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
    return 0.;
}

double StructureFunction::F2(double x, double Q2) {
    return 0.;
}

double StructureFunction::xF3(double x, double Q2) {
    return 0.;
}

double StructureFunction::F3(double x, double Q2) {
    return 0.;
}

}