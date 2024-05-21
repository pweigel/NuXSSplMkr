#include "NuXSSplMkr/structure_function.h"

namespace nuxssplmkr {

StructureFunction::StructureFunction(Configuration &_config)
    : config(_config)
{
    pc = new PhysConst();

    GF2   = SQ(pc->GF);
    M_iso = 0.5*(pc->proton_mass + pc->neutron_mass);
}

void StructureFunction::Set_Lepton_Mass(double m) {
    M_lepton = m;
}

void StructureFunction::Set_Neutrino_Energy(double E) {
    ENU = E;
}

void StructureFunction::Set_Q_APFEL(double Q) {
    if (config.SF.evolve_pdf) {
        APFEL::ComputeStructureFunctionsAPFEL(max(1.3, sqrt(config.SF.Q2min)), Q);
    } else {
        APFEL::SetAlphaQCDRef(config.pdf.pdf->alphasQ(Q), Q);
        APFEL::ComputeStructureFunctionsAPFEL(Q, Q);
    }
}

void StructureFunction::InitializeAPFEL() {

    if (config.SF.mass_scheme == "parton") {
        std::cout << "Mass scheme set to 'parton'. Not initializing APFEL!" << std::endl;
        return;
    }

    APFEL::SetPDFSet(config.pdf.pdfset);
    APFEL::SetReplica(config.pdf.replica);
    APFEL::SetMassScheme(config.SF.mass_scheme);
    std::cout << "PDF Masses: " << config.pdf.pdf_quark_masses[4] << " " << config.pdf.pdf_quark_masses[5] << " " << config.pdf.pdf_quark_masses[6] << std::endl;;
    if (config.SF.disable_top == true) { // TODO: a better way of doing this. This should only be used for CSMS I think.
        std::cout << "WARNING: Top mass set to m_b + 0.1!" << std::endl;
        APFEL::SetPoleMasses(config.pdf.pdf_quark_masses[4], config.pdf.pdf_quark_masses[5], config.pdf.pdf_quark_masses[5]+0.1);
    } else {
        APFEL::SetPoleMasses(config.pdf.pdf_quark_masses[4], config.pdf.pdf_quark_masses[5], config.pdf.pdf_quark_masses[6]);
    }

    // APFEL::SetQLimits(std::sqrt(1e-2), std::sqrt(config.SF.Q2max));
    APFEL::SetQLimits(std::sqrt(config.SF.Q2min), std::sqrt(config.SF.Q2max));
    // std::cout << "Evolution Q limits: [" << std::sqrt(config.Q2min) << ", " << std::sqrt(config.Q2max) << "] GeV." << std::endl;
    //APFEL::SetPolarizationDIS(0);
    // APFEL::EnableTargetMassCorrections(false); // Don't use this!
    std::cout << "FONLL Damping: " << config.SF.enable_FONLL_damping << ", " << config.SF.FONLL_damping_factor << std::endl;
    APFEL::EnableDampingFONLL(config.SF.enable_FONLL_damping);
    APFEL::SetDampingPowerFONLL(config.SF.FONLL_damping_factor, config.SF.FONLL_damping_factor, config.SF.FONLL_damping_factor);
    //APFEL::SetFastEvolution(true);
    //APFEL::LockGrids(true);
    //APFEL::EnableEvolutionOperator(true);
    if (config.SF.FFNS > 0) {
        APFEL::SetFFNS(config.SF.FFNS);
    }

    //APFEL::SetFFNS(3);
    //APFEL::SetTheory("QUniD");
    //APFEL::SetTheory("QED");
    //APFEL::EnableLeptonEvolution(true);
    //APFEL::SetTauMass(1e10);
    //APFEL::SetPDFEvolution("exactalpha");
    APFEL::SetNumberOfGrids(3);
    APFEL::SetGridParameters(1, 90, 3, config.SF.xmin);
    APFEL::SetGridParameters(2, 50, 5, 1e-1);
    APFEL::SetGridParameters(3, 40, 5, 8e-1);
    APFEL::SetPerturbativeOrder(config.SF.perturbative_order);
    APFEL::SetSmallxResummation(config.SF.enable_small_x, config.SF.small_x_order);
    APFEL::SetAlphaQCDRef(config.pdf.pdf->alphasQ(config.constants.MassZ), config.constants.MassZ);
    //APFEL::SetAlphaEvolution("expanded");
    //APFEL::SetPDFEvolution("expandalpha");

    APFEL::SetMaxFlavourPDFs(6);
    APFEL::SetMaxFlavourAlpha(6);
    APFEL::SetCKM(config.constants.Vud, config.constants.Vus, config.constants.Vub,
                  config.constants.Vcd, config.constants.Vcs, config.constants.Vcb,
                  config.constants.Vtd, config.constants.Vts, config.constants.Vtb);
    APFEL::SetProcessDIS(config.SF.DIS_process);

    APFEL::SetProjectileDIS(config.projectile);
    APFEL::SetTargetDIS(config.target);

    // Initializes integrals on the grids
    APFEL::InitializeAPFEL_DIS();
}

void StructureFunction::GetCoefficients() {
    int down_type[3] = {1, 3, 5};  // d s b
    int up_type[3] = {2, 4, 6};  // u c t
    // TODO: Maybe do isospin symmetry here?

    if (config.neutrino_type == neutrino) {
        for (int dt : down_type) {
            F2coef[dt] = 1.;
            F3coef[dt] = 1.;
        }
        for (int ut: up_type) {
            F2coef[-ut] = 1.;
            F3coef[-ut] = -1.;
        }
    } else if (config.neutrino_type == antineutrino) {
        for (int ut : up_type) {
            F2coef[ut] = 1.;
            F3coef[ut] = 1.;
        }
        for (int dt: down_type) {
            F2coef[-dt] = 1.;
            F3coef[-dt] = -1.;
        }
    } else {
        throw std::runtime_error("Unidentified neutrino type!");
    }

    // gluons
    F2coef[21] = 0.;
    F3coef[21] = 0.; 
}

double StructureFunction::F1(double x, double Q2) {
    if ( (config.SF.perturbative_order == LO) && (config.SF.mass_scheme == "parton") ) {
        return F2(x, Q2) / (2. * x);
    } else {
        return (F2(x, Q2) - FL(x, Q2)) / (2. * x);
    }
}

double StructureFunction::F2(double x, double Q2) {
    if ( (config.SF.mass_scheme == "parton") ) {
        auto s = PDFExtract(x, Q2);
        
        // TODO: Figure out a better way to do this (isospin symmetry for p<->n)
        if(config.target_type == neutron) {
            auto q1 = s[1];
            auto q2 = s[2];
            s[1] = q2;
            s[2] = q1;
        }
        return F2_LO(s);
    } else {
        switch(config.sf_type) {
            case SFType::total:  return APFEL::F2total(x);
            case SFType::light:  return APFEL::F2light(x);
            case SFType::charm:  return APFEL::F2charm(x);
            case SFType::bottom: return APFEL::F2bottom(x);
            case SFType::top:    return APFEL::F2top(x);
            default:             return APFEL::F2total(x);
        }
    }
}

double StructureFunction::FL(double x, double Q2) {
    if ( (config.SF.mass_scheme == "parton") ) {
        return 0.0;
    } else {
        switch(config.sf_type) {
            case SFType::total:  return APFEL::FLtotal(x);
            case SFType::light:  return APFEL::FLlight(x);
            case SFType::charm:  return APFEL::FLcharm(x);
            case SFType::bottom: return APFEL::FLbottom(x);
            case SFType::top:    return APFEL::FLtop(x);
            default:             return APFEL::FLtotal(x);
        }
    }
}

double StructureFunction::xF3(double x, double Q2) {
    if ( (config.SF.mass_scheme == "parton") ) {
        auto s = PDFExtract(x, Q2);

        // TODO: Figure out a better way to do this (isospin symmetry for p<->n)
        if(config.target_type == neutron) {
            auto q1 = s[1];
            auto q2 = s[2];
            s[1] = q2;
            s[2] = q1;
        }
        return xF3_LO(s);
    } else {
        switch(config.sf_type) {
            case SFType::total:  return APFEL::F3total(x);
            case SFType::light:  return APFEL::F3light(x);
            case SFType::charm:  return APFEL::F3charm(x);
            case SFType::bottom: return APFEL::F3bottom(x);
            case SFType::top:    return APFEL::F3top(x);
            default:             return APFEL::F3total(x);
        }
    }
}

double StructureFunction::F2_LO(map<int, double>& xq_arr) {
    /* Example: Neutrino
       F2 = 2x(d + s + b + ubar + cbar + tbar) */
    double k = 0.;
    for( int p : partons ) {
        k += F2coef[p] * xq_arr[p];
        // std::cout << p << ": " << F2coef[p] * xq_arr[p] << std::endl;
    }	
    return 2 * k;
}

double StructureFunction::xF3_LO(map<int, double>& xq_arr) {
    /* Example: Neutrino
       xF3 = 2x(d + s + b - ubar - cbar - tbar) */
    double k=0.;
    for( int p : partons ) {        
        k += F3coef[p]*xq_arr[p];
        // std::cout << k << std::endl;
    }
    
    return 2 * k;
}

double StructureFunction::F3(double x, double Q2) {
    return xF3(x, Q2) / x;
}

double StructureFunction::F4(double x, double Q2) {
    if (config.SF.use_AlbrightJarlskog) {
        return 0.0;
    } else{
        throw NotImplementedException();
    }
}

double StructureFunction::F5(double x, double Q2) {
    if (config.SF.use_AlbrightJarlskog) {
        return F2(x, Q2) / (2.0 * x);
    } else{
        throw NotImplementedException();
    }
}

double StructureFunction::RescalingVariable(double Q2) {
    double m2;
    switch(config.sf_type) {
        case SFType::total:  m2 = 0; break;
        case SFType::light:  m2 = 0; break;
        case SFType::charm:  m2 = SQ(config.pdf.pdf_quark_masses[4]); break;
        case SFType::bottom: m2 = SQ(config.pdf.pdf_quark_masses[5]); break;
        case SFType::top:    m2 = SQ(config.pdf.pdf_quark_masses[6]); break;
        default:             m2 = 0; break;
    }
    return 1.0 / (1 + m2 / Q2);
}

double StructureFunction::NachtmannR(double x, double Q2){
    double m = M_iso / (pc->GeV);
    return sqrt(1. + 4.*x*x*SQ(m)/Q2);
}

double StructureFunction::NachtmannXi(double x, double Q2){
    double denominator = 1. + NachtmannR(x, Q2);
    return 2.*x / denominator;
}

double StructureFunction::NachtmannXibar(double x, double Q2){
    double denominator = 1. + NachtmannR(x, Q2);
    double rescaling = RescalingVariable(Q2);
    return 2.*x*rescaling / denominator;
}

template<class T,double (T::*f)(double,double), int n, int m>
double StructureFunction::HGeneric(double xi, double Q2){
    // Integrate T::f from xi to 1 at Q2
    if (xi > 1.0) {
        return 0.0;
    }
    gsl_integration_cquad_workspace * w = gsl_integration_cquad_workspace_alloc(2500);    
    double result, error;
    
    gsl_function F;
    Set_Q_APFEL(std::sqrt(Q2));
    _kernel_Q2 = Q2;
    F.function = &HK<T, f, n, m>;
    F.params = this;
    size_t neval;
    // double _integrate_xmax = 1.0;
    int status = gsl_integration_cquad(&F, xi, 1.0, 0, 1.e-3, w, &result, &error, &neval);
    gsl_integration_cquad_workspace_free(w);

    return result;
}

double StructureFunction::H2(double xi, double Q2){
    return HGeneric<StructureFunction,&StructureFunction::F2,2,1>(xi, Q2);
}

double StructureFunction::H3(double xi, double Q2){
    return HGeneric<StructureFunction,&StructureFunction::F3,1,1>(xi, Q2);
}

double StructureFunction::G2(double xi, double Q2){
    return HGeneric<StructureFunction,&StructureFunction::F2,1,1>(xi, Q2) - xi*HGeneric<StructureFunction,&StructureFunction::F2,2,1>(xi, Q2);
}


double StructureFunction::F1_TMC(double x, double Q2) { // slow rescale?
    //todo: fix
    double xi = NachtmannXi(x, Q2);
    // double xibar = NachtmannXibar(x, Q2);
    if (xi > 1.0) {
        return 0.0;
    }
    double r = NachtmannR(x, Q2);
    double m = M_iso / (pc->GeV);
    double term1 = (x / (xi * r)) * F1(xi, Q2);
    double term2 = SQ(m * x / r) / Q2 * H2(xi, Q2);
    double term3 = 2.0 * SQ(SQ(m) * x / (Q2 * r)) * (x / r) * G2(xi, Q2);
    return (term1 + term2 + term3);
}

double StructureFunction::F2_TMC(double x, double Q2) {
    double xi = NachtmannXi(x, Q2);
    // double xibar = NachtmannXibar(x, Q2);
    if (xi > 1.0) {
        return 0.0;
    }
    double r = NachtmannR(x, Q2);
    double m = M_iso / (pc->GeV);

    double g2 = G2(xi, Q2);
    double h2 = H2(xi, Q2);
    
    // double term1 = SQ(x / (xi * r)) / r * F2(xi, Q2);
    // double term2 = 6.0 * SQ(m * x / SQ(r)) * (x / Q2) * H2(xi, Q2);
    // double term3 = 12.0 * SQ(SQ(m * x / r) / (Q2 * r)) * G2(xi, Q2);

    double term1 = x*x/(SQ(xi) * r*r*r) * F2(xi, Q2);
    double term2 = 6.*m*m* x*x*x  * h2 / (Q2 * r*r*r*r);
    double term3 = 12.*SQ(SQ(m*x)) * g2 / (SQ(Q2) * r*r*r*r*r );

    return (term1 + term2 + term3);
}

double StructureFunction::F3_TMC(double x, double Q2) {
    double xi = NachtmannXi(x, Q2);
    // double xibar = NachtmannXibar(x, Q2);
    if (xi > 1.0) {
      return 0.0;
    }
    double r = NachtmannR(x, Q2);
    double m = M_iso / (pc->GeV);

    // double term1 = (x / (xi * SQ(r))) * F3(xi, Q2);
    // double term2 = 2.0 * SQ(m * x / r) / (Q2 * r) * H3(xi, Q2)
    // double g2 = G2(xi, Q2);
    double h3 = H3(xi, Q2);

    double term1 = x/(xi * r*r) * F3(xi, Q2);
    double term2 = 2.*m*m* x*x  * h3 / (Q2 * r*r*r);
    
    return (term1 + term2);
}

double StructureFunction::xF3_TMC(double x, double Q2) {
    return x * F3_TMC(x, Q2);
}

double StructureFunction::CKMT_n(double Q2) {
    return (3.0 / 2.0) * (1.0 + Q2 / (Q2 + config.CKMT.c));
}

double StructureFunction::CKMT_Delta(double Q2) {
    return config.CKMT.Delta0 * (1.0 + 2.0 * Q2 / (Q2 + config.CKMT.d));
}

double StructureFunction::F_CKMT(double x, double Q2, double A, double B, double f) {
    // Jeong/Reno Paper: https://arxiv.org/pdf/2307.09241.pdf
    double delta = CKMT_Delta(Q2);
    double n = CKMT_n(Q2);
    double term1 = A * pow(x, -1.0*delta) * pow(1.0-x, n + 4.0) * pow( (Q2 / (Q2 + config.CKMT.a)), 1.0 + delta);
    double term2 = B * pow(x, 1.0 - config.CKMT.AlphaR) * pow(1.0-x, n) * pow(Q2 / (Q2 + config.CKMT.b), config.CKMT.AlphaR) * (1.0 + f * (1.0 - x));
    return (term1 + term2);
}

double StructureFunction::F2_CKMT(double x, double Q2) {
    return F_CKMT(x, Q2, config.CKMT.F2A, config.CKMT.F2B, config.CKMT.F2f);
}

double StructureFunction::xF3_CKMT(double x, double Q2) {
    return F_CKMT(x, Q2, config.CKMT.xF3A, config.cp_factor * config.CKMT.xF3B, config.CKMT.xF3f);
}

double StructureFunction::F3_CKMT(double x, double Q2) {
    return xF3_CKMT(x, Q2) / x;
}

double StructureFunction::R_CKMT(double x, double Q2){
    // R parameterization from Whitlow et al: Phys. Lett. B 250, 193 (1990)
    double big_theta = 1.0 + 12.0 * (Q2 / (Q2 + 1.0)) * (SQ(0.125) / (SQ(0.125) + SQ(x)));
    return 0.635 / log(Q2/SQ(0.2)) * big_theta + 0.5747 / Q2 - 0.3534 / (SQ(Q2) + SQ(0.3));
}

double StructureFunction::F1_CKMT(double _F2, double x, double Q2){
    double m = M_iso / pc->GeV;
    double _R = R_CKMT(x, Q2);
    // Eq. 2 from Whitlow et al:
    // R = (F2 / (2xF1)) ( 1 + 4 M^2 x^2/Q^2) - 1
    // R + 1 = 
    return _F2 * (1.0 + 4.0 * SQ(m * x) / (Q2 * SQ(pc->GeV))) / (2.0 * x * (_R + 1.0));
}

// double StructureFunction::FL_PCAC(double x, double Q2) {
//     return 0.0;
// }

double StructureFunction::F2_PCAC(double x, double Q2) {
    // Jeong/Reno Paper: https://arxiv.org/pdf/2307.09241.pdf
    double delta = CKMT_Delta(Q2);
    double n = CKMT_n(Q2);
    double M_PCAC = 0.8;
    double f_PCAC = 1.0 / SQ(1 + Q2 / SQ(M_PCAC));
    double term1 = config.PCAC.A * pow(x, -1.0*delta) * pow(1.0-x, n + 4.0) * pow( (Q2 / (Q2 + config.CKMT.a)), delta);
    double term2 = config.PCAC.B * pow(x, 1.0 - config.CKMT.AlphaR) * pow(1.0-x, n) * pow(Q2 / (Q2 + config.CKMT.b), config.CKMT.AlphaR - 1.0) * (1.0 + config.CKMT.F2f * (1.0 - x));
    return f_PCAC * (term1 + term2);
}

std::map<int,double> StructureFunction::PDFExtract(double x, double Q2){
    LHAPDF::GridPDF* grid_central = dynamic_cast<LHAPDF::GridPDF*>(config.pdf.pdf);
    string xt = "nearest";
    grid_central->setExtrapolator(xt);

    std::map<int,double> xq_arr;
    for ( int p : partons ){
      xq_arr[p] = grid_central -> xfxQ2(p, x, Q2/(pc->GeV2));
    }
    return xq_arr;
}

std::tuple<double,double,double> StructureFunction::EvaluateSFs(double x, double Q2) {
    double _FL, _F1, _F2, _F3;

    double Q2_eval = Q2;  // Q2 that the SFs are evaluated at
    if ( (config.SF.enable_CKMT) ) {
        Q2_eval = SQ(config.CKMT.Q0); // get CKMT reference Q0
    }

    if ((config.SF.mass_scheme != "parton") || (Q2_eval != _Q2_cached)) {
        Set_Q_APFEL(std::sqrt(Q2_eval));
        _Q2_cached = Q2_eval;  // We cache this so we don't keep running the APFEL function
    }

    // Get F1/F2/F3
    if ( (config.SF.enable_TMC ) && (Q2 < config.SF.TMC_Q2max)) {
        _F1 = F1_TMC(x, Q2_eval);
        _F2 = F2_TMC(x, Q2_eval);
        _F3 = F3_TMC(x, Q2_eval);
    } else {
        _FL = FL(x, Q2_eval);
        _F2 = F2(x, Q2_eval);
        _F1 = (_F2 - _FL) / (2. * x);
        _F3 = F3(x, Q2_eval);
    }

    if (config.SF.enable_CKMT) {
        /* 
        When using CKMT, we evaluate the APFEL-based SFs using the threshold
        value of Q0 and use the parameterized SFs below the threshold.
        */
        double CKMT_Q2 = SQ(config.CKMT.Q0);
        // std::cout << "Q2_eval = " << Q2_eval << ", Q2 = " << Q2 << ", CKMT_Q2 = " << CKMT_Q2 << std::endl;

        if (Q2 < CKMT_Q2) {
            double _F2_CKMT    = F2_CKMT(x, Q2);
            double _F3_CKMT    = F3_CKMT(x, Q2);
            double _F1_CKMT    = F1_CKMT(_F2_CKMT, x, Q2);
            
            double _F1_Q0      = _F1;
            double _F2_Q0      = _F2;
            double _F3_Q0      = _F3;
            double _F2_CKMT_Q0 = F2_CKMT(x, Q2_eval);
            double _F3_CKMT_Q0 = F3_CKMT(x, Q2_eval);
            double _F1_CKMT_Q0 = F1_CKMT(_F2_CKMT_Q0, x, Q2_eval);

            double _F2_PCAC;
            double _F2_PCAC_Q0;
            if (config.SF.enable_PCAC) {
                _F2_PCAC = F2_PCAC(x, Q2);
                _F2_PCAC_Q0 = F2_PCAC(x, Q2_eval);
                _F2_CKMT += _F2_PCAC;
                _F2_CKMT_Q0 += _F2_PCAC_Q0;
            }

            _F2 = _F2_CKMT * (_F2_Q0 / _F2_CKMT_Q0);
            _F1 = _F1_CKMT * (_F1_Q0 / _F1_CKMT_Q0);
            _F3 = _F3_CKMT * (_F3_Q0 / _F3_CKMT_Q0);
        } else {
            // if Q^2 > Q0^2, then we use the regular SFs (do nothing)
        }
    }

    return {_F1, _F2, _F3};
}

void StructureFunction::Compute() {
    /*
    Refactored version of the calculations fron BuildSplines
    */
    const unsigned int Nx = config.SF.Nx;
    const unsigned int NQ2 = config.SF.NQ2;

    std::vector<double> x_arr;
    std::vector<double> Q2_arr;

    // Get the coefficients for parton calculation
    if (config.SF.mass_scheme == "parton") {
        GetCoefficients();
    }

    // Step sizes in log space
    // TODO: log-lin spacing in x at some point
    double d_log_Q2 = std::abs( std::log10(config.SF.Q2min) - std::log10(config.SF.Q2max) ) / (NQ2 - 1);
    double d_log_x  = std::abs( std::log10(config.SF.xmin)  - std::log10(config.SF.xmax)  ) / (Nx - 1);

    // Get the Q2 and x values
    size_t N_samples = Nx * NQ2;
    for ( double log_Q2 = std::log10(config.SF.Q2min); log_Q2<std::log10(config.SF.Q2max); log_Q2 += d_log_Q2 ) {
        double powQ2 = std::pow( 10, log_Q2);
        double Q2 = log_Q2;
        if (powQ2 > config.SF.Q2max) continue;
        Q2_arr.push_back(Q2);
    }

    for ( double log_x = std::log10(config.SF.xmin); log_x<std::log10(config.SF.xmax); log_x += d_log_x ) {
        double powx = std::pow( 10, log_x);
        double x = log_x ;
        if ( powx > config.SF.xmax ) continue;
        x_arr.push_back(x);
    }

    std::tuple<double,double,double> sfs;
    for (unsigned int Q2i = 0; Q2i < NQ2; Q2i++) {
        double log_Q2 = std::log10(config.SF.Q2min) + Q2i * d_log_Q2;
        double Q2 = std::pow(10.0, log_Q2);

        for (unsigned int xi = 0; xi < Nx; xi++) {
            double log_x = std::log10(config.SF.xmin) + xi * d_log_x;
            double x = std::pow(10.0, log_x);
            sfs = EvaluateSFs(x, Q2);
        }
    }
}

void StructureFunction::LoadGrids(string inpath) {

}

void StructureFunction::BuildGrids(string outpath) {
    save_splines = false;
    BuildSplines(outpath);
}

// TODO: The order of x, Q2 in splines is not consistent with the function definitions here
void StructureFunction::BuildSplines(string outpath) {
    const unsigned int Nx = config.SF.Nx;
    const unsigned int NQ2 = config.SF.NQ2;

    std::vector<double> x_arr;
    std::vector<double> Q2_arr;

    // Get the coefficients for parton calculation
    if (config.SF.mass_scheme == "parton") {
        GetCoefficients();
    }

    // Step sizes in log space
    double d_log_Q2 = std::abs( std::log10(config.SF.Q2min) - std::log10(config.SF.Q2max) ) / (NQ2 - 1);
    double d_log_x  = std::abs( std::log10(config.SF.xmin)  - std::log10(config.SF.xmax)  ) / (Nx - 1);

    // Spline parameters
    const uint32_t dim = 2;
    std::vector<uint32_t> orders(dim, 2);

    unsigned int Nknots_Q2 = config.SF.Nx + 2;
    unsigned int Nknots_x = config.SF.NQ2 + 2;

    std::vector<double> Q2_knots;
    std::vector<double> x_knots;

    // Step sizes for knots in log space
    const double d_log_Q2_knot= std::abs(std::log10(config.SF.Q2max)- std::log10(config.SF.Q2min)) / (Nknots_Q2 - 1);
    const double d_log_x_knot = std::abs(std::log10(config.SF.xmax) - std::log10(config.SF.xmin) ) / (Nknots_x - 1);

    std::cout << "log_Q2min = " << std::log10(config.SF.Q2min) << ", log_Q2max = " << std::log10(config.SF.Q2max) << std::endl;
    std::cout << "log_xmin = " << std::log10(config.SF.xmin) << ", log_xmax = " << std::log10(config.SF.xmax) << std::endl;
    std::cout << "d_log_Q2 = " << d_log_Q2 << ", d_log_x = " << d_log_x << std::endl;

    // Get the Q2 and x values
    size_t N_samples = Nx * NQ2;
    for ( double log_Q2 = std::log10(config.SF.Q2min); log_Q2<std::log10(config.SF.Q2max); log_Q2 += d_log_Q2 ) {
        double powQ2 = std::pow( 10, log_Q2);
        double Q2 = log_Q2;
        if (powQ2 > config.SF.Q2max) continue;
        Q2_arr.push_back(Q2);
    }

    for ( double log_x = std::log10(config.SF.xmin); log_x<std::log10(config.SF.xmax); log_x += d_log_x ) {
        double powx = std::pow( 10, log_x);
        double x = log_x ;
        if ( powx > config.SF.xmax ) continue;
        x_arr.push_back(x);
    }

    // Get the knots for Q2 and x
    // for ( double log_Q2 = std::log10(config.SF.Q2min) - d_log_Q2_knot; log_Q2 <= std::log10(config.SF.Q2max) + d_log_Q2_knot; log_Q2 += d_log_Q2_knot ) {
    //     double knot = log_Q2;
    //     Q2_knots.push_back(knot);
    // }
    // for ( double log_x = std::log10(config.SF.xmin) - d_log_x_knot; log_x <= 1 + d_log_x_knot; log_x += d_log_x_knot ) {
    //     double knot = log_x;
    //     x_knots.push_back(knot);
    // }

    // ### Testing new knots ###
    for ( double log_Q2 = 0; log_Q2 <= 1; log_Q2 += 0.01 ) {
        double knot = log_Q2;
        Q2_knots.push_back(knot);
    }

    for ( double log_Q2 = 1 + d_log_Q2_knot; log_Q2 <= std::log10(config.SF.Q2max) + d_log_Q2_knot; log_Q2 += d_log_Q2_knot ) {
        double knot = log_Q2;
        Q2_knots.push_back(knot);
    }

    for ( double log_x = std::log10(config.SF.xmin) - d_log_x_knot; log_x <= -1.0-d_log_x_knot; log_x += d_log_x_knot ) {
        double knot = log_x;
        x_knots.push_back(knot);
    }

    for ( double log_x = -1.0; log_x <= 0.25; log_x += 0.01 ) {
        double knot = log_x;
        x_knots.push_back(knot);
    }
    // ####

    // setup grid stuff
    string f1_grid_fn = outpath + "/F1_"+config.projectile+"_"+config.target+"_"+config.sf_type_string+".grid";
    string f2_grid_fn = outpath + "/F2_"+config.projectile+"_"+config.target+"_"+config.sf_type_string+".grid";
    string f3_grid_fn = outpath + "/F3_"+config.projectile+"_"+config.target+"_"+config.sf_type_string+".grid";

    std::ofstream f1_outfile;
    f1_outfile.open(f1_grid_fn);
    std::ofstream f2_outfile;
    f2_outfile.open(f2_grid_fn);
    std::ofstream f3_outfile;
    f3_outfile.open(f3_grid_fn);

    // Write header
    f1_outfile << NQ2 << " " << Nx << "\n";
    f1_outfile << std::log10(config.SF.Q2min) << " " << std::log10(config.SF.Q2max) << " " << std::log10(config.SF.xmin) << " " << std::log10(config.SF.xmax) << "\n";
    f2_outfile << NQ2 << " " << Nx << "\n";
    f2_outfile << std::log10(config.SF.Q2min) << " " << std::log10(config.SF.Q2max) << " " << std::log10(config.SF.xmin) << " " << std::log10(config.SF.xmax) << "\n";
    f3_outfile << NQ2 << " " << Nx << "\n";
    f3_outfile << std::log10(config.SF.Q2min) << " " << std::log10(config.SF.Q2max) << " " << std::log10(config.SF.xmin) << " " << std::log10(config.SF.xmax) << "\n";
    

    // Collect SF values and write grids
    std::deque<std::pair<double,std::array<unsigned int, 2>>> F1_spline_data;
    std::deque<std::pair<double,std::array<unsigned int, 2>>> F2_spline_data;
    std::deque<std::pair<double,std::array<unsigned int, 2>>> F3_spline_data;
    for (unsigned int Q2i = 0; Q2i < NQ2; Q2i++) {
        double log_Q2 = std::log10(config.SF.Q2min) + Q2i * d_log_Q2;
        double Q2 = std::pow(10.0, log_Q2);
        if (config.general.debug) {
            std::cout << "Q2 = " << Q2 << std::endl;
        }

        double Q2eval = Q2;
        if ( (config.SF.enable_CKMT) && (Q2 < SQ(config.CKMT.Q0))) {
            Q2eval = SQ(config.CKMT.Q0);
        }

        if (config.SF.mass_scheme != "parton") {
            Set_Q_APFEL(std::sqrt(Q2eval));
        }

        for (unsigned int xi = 0; xi < Nx; xi++) {

            double log_x = std::log10(config.SF.xmin) + xi * d_log_x;
            double x = std::pow(10.0, log_x);

            double _FL, _F1, _F2, _F3;
            
            std::tie(_F1, _F2, _F3) = EvaluateSFs(x, Q2);

            if(!std::isfinite(_F1)) {
                std::cerr << "F1 Infinite! Q2 = " << Q2 << ", x = " << x << ". Setting to zero." << std::endl;
                _F1 = 0.0;
            }
            if(!std::isfinite(_F2)) {
                std::cerr << "F2 Infinite! Q2 = " << Q2 << ", x = " << x << ". Setting to zero." << std::endl;
                _F2 = 0.0;
            }
            if(!std::isfinite(_F3)) {
                std::cerr << "F3 Infinite! Q2 = " << Q2 << ", x = " << x << ". Setting to zero." << std::endl;
                _F3 = 0.0;
            }

            // Collect spline data
            F1_spline_data.push_back(std::make_pair(_F1, std::array<unsigned int, 2>{Q2i, xi}));
            F2_spline_data.push_back(std::make_pair(_F2, std::array<unsigned int, 2>{Q2i, xi}));
            F3_spline_data.push_back(std::make_pair(_F3, std::array<unsigned int, 2>{Q2i, xi}));

            // Write grid data
            f1_outfile << _F1;
            f2_outfile << _F2;
            f3_outfile << _F3;
            if (xi < Nx-1) {
              f1_outfile << ",";
              f2_outfile << ",";
              f3_outfile << ",";
            } else {
              f1_outfile << "\n";
              f2_outfile << "\n";
              f3_outfile << "\n";
            }
        }
    }

    if (save_splines) {
        // Put data into ndsparses
        photospline::ndsparse F1_data(F1_spline_data.size(), 2);
        for(auto& entry : F1_spline_data) {
            F1_data.insertEntry(entry.first, &entry.second[0]);
        }
        photospline::ndsparse F2_data(F2_spline_data.size(), 2);
        for(auto& entry : F2_spline_data) {
            F2_data.insertEntry(entry.first, &entry.second[0]);
        }
        photospline::ndsparse F3_data(F3_spline_data.size(), 2);
        for(auto& entry : F3_spline_data) {
            F3_data.insertEntry(entry.first, &entry.second[0]);
        }

        // TODO: These should all be the same size
        std::vector<double> F1_weights(F1_spline_data.size(),1.);
        std::vector<double> F2_weights(F2_spline_data.size(),1.);
        std::vector<double> F3_weights(F3_spline_data.size(),1.);

        // double smooth_Q2 = 1;
        // double smooth_x = 1;
        double smooth = 1e-5;

        // Fit splines
        photospline::splinetable<> F1_spline;
        F1_spline.fit(F1_data, F1_weights, std::vector<std::vector<double>>{Q2_arr, x_arr}, orders, 
                      {Q2_knots, x_knots}, {smooth, smooth}, {2, 2});
        if(std::isnan(*F1_spline.get_coefficients())){
            std::cerr << "F1 spline fit has failed!" << std::endl;
        }

        photospline::splinetable<> F2_spline;
        F2_spline.fit(F2_data, F2_weights, std::vector<std::vector<double>>{Q2_arr, x_arr}, orders, 
                      {Q2_knots, x_knots}, {smooth, smooth}, {2, 2});
        if(std::isnan(*F2_spline.get_coefficients())){
            std::cerr << "F2 spline fit has failed!" << std::endl;
        }

        photospline::splinetable<> F3_spline;
        F3_spline.fit(F3_data, F3_weights, std::vector<std::vector<double>>{Q2_arr, x_arr}, orders, 
                      {Q2_knots, x_knots}, {smooth, smooth}, {2, 2});
        if(std::isnan(*F3_spline.get_coefficients())){
            std::cerr << "F3 spline fit has failed!" << std::endl;
        }

        // Write splines
        F1_spline.write_fits(outpath + "/F1_"+config.projectile+"_"+config.target+"_"+config.sf_type_string+".fits");
        F2_spline.write_fits(outpath + "/F2_"+config.projectile+"_"+config.target+"_"+config.sf_type_string+".fits");
        F3_spline.write_fits(outpath + "/F3_"+config.projectile+"_"+config.target+"_"+config.sf_type_string+".fits");
    }
}


void StructureFunction::BuildGrids_v2(string outpath) {
    const unsigned int Nx = config.SF.Nx;
    const unsigned int NQ2 = config.SF.NQ2;

    std::vector<double> x_arr;
    std::vector<double> Q2_arr;

    // Get the coefficients for parton calculation
    if (config.SF.mass_scheme == "parton") {
        GetCoefficients();
    }

    // Step sizes in log space
    double d_log_Q2 = std::abs( std::log10(config.SF.Q2min) - std::log10(config.SF.Q2max) ) / (NQ2-1);
    double d_log_x  = std::abs( std::log10(config.SF.xmin)  - std::log10(config.SF.xmax)  ) / (Nx-1);

    std::cout << "log_Q2min = " << std::log10(config.SF.Q2min) << ", log_Q2max = " << std::log10(config.SF.Q2max) << std::endl;
    std::cout << "log_xmin = " << std::log10(config.SF.xmin) << ", log_xmax = " << std::log10(config.SF.xmax) << std::endl;
    std::cout << "d_log_Q2 = " << d_log_Q2 << ", d_log_x = " << d_log_x << std::endl;

    // setup grid stuff
    string f1_grid_fn = outpath + "/F1_"+config.projectile+"_"+config.target+"_"+config.sf_type_string+".grid";
    string f2_grid_fn = outpath + "/F2_"+config.projectile+"_"+config.target+"_"+config.sf_type_string+".grid";
    string f3_grid_fn = outpath + "/F3_"+config.projectile+"_"+config.target+"_"+config.sf_type_string+".grid";

    std::ofstream f1_outfile;
    f1_outfile.open(f1_grid_fn);
    std::ofstream f2_outfile;
    f2_outfile.open(f2_grid_fn);
    std::ofstream f3_outfile;
    f3_outfile.open(f3_grid_fn);

    // Get the Q2 and x values
    size_t N_samples = Nx * NQ2;
    for (int i = 0; i < NQ2; i++) {
        double log_Q2 = std::log10(config.SF.Q2min) + i * d_log_Q2;
        Q2_arr.push_back(log_Q2);
        f1_outfile << log_Q2 << " "; f2_outfile << log_Q2 << " "; f3_outfile << log_Q2 << " ";
    }
    f1_outfile << "\n"; f2_outfile << "\n"; f3_outfile << "\n";

    for (int j = 0; j < Nx; j++) {
        double log_x = std::log10(config.SF.xmin) + j * d_log_x;
        x_arr.push_back(log_x);
        f1_outfile << log_x << " "; f2_outfile << log_x << " "; f3_outfile << log_x << " ";
    }
    f1_outfile << "\n"; f2_outfile << "\n"; f3_outfile << "\n";

    // for ( double log_Q2 = std::log10(config.SF.Q2min); log_Q2<=std::log10(config.SF.Q2max); log_Q2 += d_log_Q2 ) {
    //     double powQ2 = std::pow( 10, log_Q2);
    //     double Q2 = log_Q2;
    //     Q2_arr.push_back(Q2);
    //     f1_outfile << Q2 << " "; f2_outfile << Q2 << " "; f3_outfile << Q2 << " ";
    // }

    // for ( double log_x = std::log10(config.SF.xmin); log_x<=std::log10(config.SF.xmax); log_x += d_log_x ) {
    //     double powx = std::pow( 10, log_x);
    //     double x = log_x ;
    //     x_arr.push_back(x);
    //     f1_outfile << x << " "; f2_outfile << x << " "; f3_outfile << x << " ";
    // }
    // f1_outfile << "\n"; f2_outfile << "\n"; f3_outfile << "\n";

    // Collect SF values and write grids
    for (unsigned int i = 0; i < NQ2; i++) {
        double log_Q2 = Q2_arr.at(i);
        double Q2 = std::pow(10.0, log_Q2);
        if (config.general.debug) {
            std::cout << "Q2 = " << Q2 << std::endl;
        }

        double Q2eval = Q2;
        if ( (config.SF.enable_CKMT) && (Q2 < SQ(config.CKMT.Q0))) {
            Q2eval = SQ(config.CKMT.Q0);
        }

        if (config.SF.mass_scheme != "parton") {
            Set_Q_APFEL(std::sqrt(Q2eval));
        }

        for (unsigned int j = 0; j < Nx; j++) {
            double log_x = x_arr.at(j);
            double x = std::pow(10.0, log_x);

            double _FL, _F1, _F2, _F3;
            std::tie(_F1, _F2, _F3) = EvaluateSFs(x, Q2);

            if(!std::isfinite(_F1)) {
                std::cerr << "F1 Infinite! Q2 = " << Q2 << ", x = " << x << ". Setting to zero." << std::endl;
                _F1 = 0.0;
            }
            if(!std::isfinite(_F2)) {
                std::cerr << "F2 Infinite! Q2 = " << Q2 << ", x = " << x << ". Setting to zero." << std::endl;
                _F2 = 0.0;
            }
            if(!std::isfinite(_F3)) {
                std::cerr << "F3 Infinite! Q2 = " << Q2 << ", x = " << x << ". Setting to zero." << std::endl;
                _F3 = 0.0;
            }

            // Write grid data
            f1_outfile << _F1;
            f2_outfile << _F2;
            f3_outfile << _F3;
            if (j < Nx-1) {
              f1_outfile << ",";
              f2_outfile << ",";
              f3_outfile << ",";
            } else {
              f1_outfile << "\n";
              f2_outfile << "\n";
              f3_outfile << "\n";
            }
        }
    }
}

}