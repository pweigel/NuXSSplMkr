#include "NuXSSplMkr/structure_function.h"

namespace nuxssplmkr {

StructureFunction::StructureFunction(Configuration &_config)
    : config(_config)
{
    pc = new PhysConst();

    GF2   = SQ(pc->GF);

    F1_code = std::stoi(config.Get_SF_Code("F1"));
    F2_code = std::stoi(config.Get_SF_Code("F2"));
    F3_code = std::stoi(config.Get_SF_Code("F3"));

    Set_Mode(config.mode);
}

void StructureFunction::Set_Mode(int _mode) {
    mode = _mode;

    switch (mode) {
        case 0: {insuffix = ""; outsuffix = ""; break;}
        case 1: {insuffix = ""; outsuffix = ""; break;}
        case 2: {insuffix = ""; outsuffix = "_TMC"; break;}
        case 3: {insuffix = "_TMC"; outsuffix = "_CKMT"; break;}
        case 4: {insuffix = "_TMC"; outsuffix = "_PCAC"; break;}
    }
}

void StructureFunction::Set_Neutrino_Energy(double E) {
    ENU = E;
}

void StructureFunction::Set_Q_APFEL(double Q) {
    if (config.SF.evolve_pdf) {
        APFEL::ComputeStructureFunctionsAPFEL(max(1.3, sqrt(config.SF.Q2min)), Q);
    } else {
        APFEL::SetAlphaQCDRef(config.pdf->alphasQ(Q), Q);
        APFEL::ComputeStructureFunctionsAPFEL(Q, Q);
    }
}

void StructureFunction::InitializeAPFEL() {

    if (config.SF.mass_scheme == "parton") {
        std::cout << "Mass scheme set to 'parton'. Not initializing APFEL!" << std::endl;
        return;
    }

    APFEL::SetPDFSet(config.pdf_info.pdfset);
    APFEL::SetReplica(config.pdf_info.replica);
    APFEL::SetMassScheme(config.SF.mass_scheme);
    std::cout << "PDF Masses: " << config.pdf_info.pdf_quark_masses[4] << " " << config.pdf_info.pdf_quark_masses[5] << " " << config.pdf_info.pdf_quark_masses[6] << std::endl;;
    if (config.SF.disable_top == true) { // TODO: a better way of doing this. This should only be used for CSMS I think.
        // std::cout << "WARNING: Top mass set to m_b + 0.1!" << std::endl;
        std::cout << "WARNING: Top mass changed to disable it!" << std::endl;
        APFEL::SetPoleMasses(config.pdf_info.pdf_quark_masses[4], config.pdf_info.pdf_quark_masses[5], 1e12);
    } else {
        APFEL::SetPoleMasses(config.pdf_info.pdf_quark_masses[4], config.pdf_info.pdf_quark_masses[5], config.pdf_info.pdf_quark_masses[6]);
    }

    APFEL::SetQLimits(std::sqrt(config.SF.Q2min), std::sqrt(config.SF.Q2max));
    
    std::cout << "FONLL Damping: " << config.SF.enable_FONLL_damping << ", " << config.SF.FONLL_damping_factor << std::endl;
    APFEL::EnableDampingFONLL(config.SF.enable_FONLL_damping);
    APFEL::SetDampingPowerFONLL(config.SF.FONLL_damping_factor, config.SF.FONLL_damping_factor, config.SF.FONLL_damping_factor);

    if (config.SF.FFNS > 0) {
        APFEL::SetFFNS(config.SF.FFNS);
        // APFEL::SetVFNS();
    }

    APFEL::SetNumberOfGrids(1);
    APFEL::SetGridParameters(1, 190, 3, config.SF.xmin);
    // APFEL::SetGridParameters(2, 40, 5, 1e-1);
    // APFEL::SetGridParameters(3, 20, 5, 8e-1);
    APFEL::SetPerturbativeOrder(config.SF.perturbative_order);
    APFEL::SetSmallxResummation(config.SF.enable_small_x, config.SF.small_x_order);
    APFEL::SetAlphaQCDRef(config.pdf->alphasQ(config.constants.MassZ), config.constants.MassZ);

    APFEL::SetMaxFlavourPDFs(config.SF.nf);
    APFEL::SetMaxFlavourAlpha(config.SF.nf);
    APFEL::SetCKM(config.constants.Vud, config.constants.Vus, config.constants.Vub,
                  config.constants.Vcd, config.constants.Vcs, config.constants.Vcb,
                  config.constants.Vtd, config.constants.Vts, config.constants.Vtb);

    APFEL::SetProcessDIS(config.current);

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
        double m = (config.target_mass) / (pc->GeV);
        return (F2(x, Q2) * (1.0 + 4.0*SQ(x*m)/Q2)  - FL(x, Q2)) / (2. * x);
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

double StructureFunction::RescalingVariable(double x, double Q2) {
    double m2;
    switch(config.sf_type) {
        case SFType::total:  m2 = 0; break;
        case SFType::light:  m2 = 0; break;
        case SFType::charm:  m2 = SQ(config.pdf_info.pdf_quark_masses[4]); break;
        case SFType::bottom: m2 = SQ(config.pdf_info.pdf_quark_masses[5]); break;
        case SFType::top:    m2 = SQ(config.pdf_info.pdf_quark_masses[6]); break;
        default:             m2 = 0; break;
    }
    return x * (1 + m2 / Q2);
}

double StructureFunction::NachtmannR(double x, double Q2){
    double m = (config.target_mass) / (pc->GeV);
    return sqrt(1. + 4.*x*x*SQ(m)/Q2);
}

double StructureFunction::NachtmannXi(double x, double Q2){
    double denominator = 1. + NachtmannR(x, Q2);
    return 2.*x / denominator;
}

// double StructureFunction::NachtmannXibar(double x, double Q2){
//     double denominator = 1. + NachtmannR(x, Q2);
//     double rescaling = RescalingVariable(Q2);
//     return 2.*x*rescaling / denominator;
// }

template<class T,double (T::*f)(double,double), int n, int m>
double StructureFunction::HGeneric(double xi, double Q2){
    // Integrate T::f from xi to 1 at Q2
    if ( (xi > 1.0) || (xi < config.SF.xmin)) {
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
    gsl_integration_cquad(&F, xi, 1.0, 0, 1.e-3, w, &result, &error, &neval);
    gsl_integration_cquad_workspace_free(w);

    return result;
}

double StructureFunction::H2_kernel(double u) {
    double f2 = SF_PDF->xfxQ2(F2_code, u, _kernel_Q2);
    return f2 / SQ(u);
}

double StructureFunction::H3_kernel(double u) {
    double f3 = SF_PDF->xfxQ2(F2_code, u, _kernel_Q2);
    return f3 / u;
}

double StructureFunction::G2_kernel(double v) {
    double f2 = SF_PDF->xfxQ2(F2_code, v, _kernel_Q2);
    return (v - _kernel_xi) * f2 / SQ(v);
}

double StructureFunction::H2(double xi, double Q2) {
    if ( (xi > 1.0) || (xi < config.SF.xmin)) {
        return 0.0;
    }

    gsl_integration_cquad_workspace * w = gsl_integration_cquad_workspace_alloc(2500);    
    double result, error;
    
    gsl_function F;
    _kernel_Q2 = Q2;

    F.function = &KernelWrapper<StructureFunction, &StructureFunction::H2_kernel>;
    F.params = this;
    size_t neval;

    gsl_integration_cquad(&F, xi, 1.0, 0, 1.e-5, w, &result, &error, &neval);
    gsl_integration_cquad_workspace_free(w);

    return result;
}

double StructureFunction::H3(double xi, double Q2) {
    if ( (xi > 1.0) || (xi < config.SF.xmin)) {
        return 0.0;
    }

    gsl_integration_cquad_workspace * w = gsl_integration_cquad_workspace_alloc(2500);    
    double result, error;
    
    gsl_function F;
    _kernel_Q2 = Q2;

    F.function = &KernelWrapper<StructureFunction, &StructureFunction::H3_kernel>;
    F.params = this;
    size_t neval;

    gsl_integration_cquad(&F, xi, 1.0, 0, 1.e-5, w, &result, &error, &neval);
    gsl_integration_cquad_workspace_free(w);

    return result;
}

double StructureFunction::G2(double xi, double Q2) {
    if ( (xi > 1.0) || (xi < config.SF.xmin)) {
        return 0.0;
    }

    gsl_integration_cquad_workspace * w = gsl_integration_cquad_workspace_alloc(2500);    
    double result, error;
    
    gsl_function F;
    _kernel_Q2 = Q2;
    _kernel_xi = xi;

    F.function = &KernelWrapper<StructureFunction, &StructureFunction::G2_kernel>;
    F.params = this;
    size_t neval;

    gsl_integration_cquad(&F, xi, 1.0, 0, 1.e-5, w, &result, &error, &neval);
    gsl_integration_cquad_workspace_free(w);

    return result;
}

double StructureFunction::F1_TMC(double x, double Q2) {
    double xi = NachtmannXi(x, Q2);
    if ( (xi > 1.0) || (xi < config.SF.xmin)) {
        return 0.0;
    }
    double r = NachtmannR(x, Q2);
    double m = (config.target_mass) / (pc->GeV);

    // Get the value from the F1 spline at Q2, xi
    double f1 = SF_PDF->xfxQ2(F1_code, xi, Q2);

    double h2 = H2(xi, Q2);
    double g2 = G2(xi, Q2);

    double term1 = (x / (xi * r)) * f1;
    double term2 = SQ(m * x / r) / Q2 * h2;
    double term3 = 2.0 * SQ(SQ(m) * x / (Q2 * r)) * (x / r) * g2;

    return (term1 + term2 + term3);
}

double StructureFunction::F2_TMC(double x, double Q2) {
    double xi = NachtmannXi(x, Q2);
    if ( (xi > 1.0) || (xi < config.SF.xmin)) {
        return 0.0;
    }

    double r = NachtmannR(x, Q2);
    double m = (config.target_mass) / (pc->GeV);

    // Get the value from the F2 spline at Q2, xi
    double f2 = SF_PDF->xfxQ2(F2_code, xi, Q2);

    double g2 = G2(xi, Q2);
    double h2 = H2(xi, Q2);

    double term1 = x*x/(SQ(xi) * r*r*r) * f2;
    double term2 = 6.*m*m* x*x*x  * h2 / (Q2 * r*r*r*r);
    double term3 = 12.*SQ(SQ(m*x)) * g2 / (SQ(Q2) * r*r*r*r*r );

    return (term1 + term2 + term3);
}

double StructureFunction::F3_TMC(double x, double Q2) {
    double xi = NachtmannXi(x, Q2);
    // double xibar = NachtmannXibar(x, Q2);
    if ( (xi > 1.0) || (xi < config.SF.xmin)) {
      return 0.0;
    }
    double r = NachtmannR(x, Q2);
    double m = (config.target_mass) / (pc->GeV);

    // Get the value from the F2 spline at Q2, xi
    double f3 = SF_PDF->xfxQ2(F3_code, xi, Q2);

    double h3 = H3(xi, Q2);

    double term1 = x/(xi * r*r) * f3;
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
    // Note: The value of alphaR in this paper may be a type, the old paper (hep-ph 0605295) uses 0.4250 instead of 0.4150
    double delta = CKMT_Delta(Q2);
    double n = CKMT_n(Q2);
    double term1 = A * pow(x, -1.0*delta) * pow(1.0-x, n + 4.0) * pow( (Q2 / (Q2 + config.CKMT.a)), 1.0 + delta);
    double term2 = B * pow(x, 1.0 - config.CKMT.AlphaR) * pow(1.0-x, n) * pow(Q2 / (Q2 + config.CKMT.b), config.CKMT.AlphaR) * (1.0 + f * (1.0 - x));
    return (term1 + term2);
}

double StructureFunction::F2_CKMT(double x, double Q2) {
    if (config.sf_type == SFType::charm) {
        return F_CKMT(x, Q2, config.CKMT.F2A, 0.0, config.CKMT.F2f); // sea only
    }
    return F_CKMT(x, Q2, config.CKMT.F2A, config.CKMT.F2B, config.CKMT.F2f);
}

double StructureFunction::xF3_CKMT(double x, double Q2) {
    if (config.sf_type == SFType::charm) {
        return F_CKMT(x, Q2, config.CKMT.xF3A, 0.0, config.CKMT.xF3f);
    }
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
    double m = (config.target_mass) / (pc->GeV);
    double _R;
    if (Q2 < 0.3) {
        // _R = R_CKMT(x, 0.3) * Q2 / 0.3;
        _R = R_CKMT(x, 0.3) * 0.3 / Q2;
    } else {
        _R = R_CKMT(x, Q2);
    }
    // Eq. 2 from Whitlow et al:
    // R = (F2 / (2xF1)) ( 1 + 4 M^2 x^2/Q^2) - 1
    // return _F2 * (1.0 + 4.0 * SQ(m * x) / (Q2 * SQ(pc->GeV))) / (2.0 * x * (_R + 1.0));
    return _F2 * (1.0 + 4.0 * SQ(m * x) / (Q2)) / (2.0 * x * (_R + 1.0));
}

double StructureFunction::F2_PCAC(double x, double Q2) {
    // Jeong/Reno Paper: https://arxiv.org/pdf/2307.09241.pdf
    double delta = CKMT_Delta(Q2);
    double n = CKMT_n(Q2);
    double M_PCAC = 0.8; // GeV
    double f_PCAC = 1.0 / SQ(1 + Q2 / SQ(M_PCAC));
    double term1 = config.PCAC.A * pow(x, -1.0*delta) * pow(1.0-x, n + 4.0) * pow( (Q2 / (Q2 + config.CKMT.a)), delta);
    double term2 = config.PCAC.B * pow(x, 1.0 - config.CKMT.AlphaR) * pow(1.0-x, n) * pow(Q2 / (Q2 + config.CKMT.b), config.CKMT.AlphaR - 1.0) * (1.0 + config.CKMT.F2f * (1.0 - x));
    if (config.sf_type == SFType::charm) {
        term2 = 0.;
    }
    return f_PCAC * (term1 + term2);
}

std::map<int,double> StructureFunction::PDFExtract(double x, double Q2){
    LHAPDF::GridPDF* grid_central = dynamic_cast<LHAPDF::GridPDF*>(config.pdf);
    string xt = "nearest";
    grid_central->setExtrapolator(xt);

    std::map<int,double> xq_arr;
    for ( int p : partons ){
      xq_arr[p] = grid_central -> xfxQ2(p, x, Q2/(pc->GeV2));
    }
    return xq_arr;
}

std::tuple<double,double,double,double> StructureFunction::EvaluateSFs(double x, double Q2) {
    double _F1 = 0, _F2 = 0, _F3 = 0, _FL = 0;

    switch (mode) {
        case 0:  // Parton
            break;
        case 1: // Plain ol' APFEL structure functions
        {
            double m = (config.target_mass) / (pc->GeV);
            _FL = FL(x, Q2);
            _F2 = F2(x, Q2);
            // _F1 = (_F2 * (1.0 + 4.0*SQ(x*m)/Q2) - _FL) / (2. * x);
            _F1 = (_F2 - _FL) / (2. * x);
            _F3 = F3(x, Q2);
            break;
        }
        case 2: // TMC
        {
            if (x > 1e-3) {
                _F1 = F1_TMC(x, Q2);
                _F2 = F2_TMC(x, Q2);
                _F3 = F3_TMC(x, Q2);
            } else {
                _F1 = SF_PDF->xfxQ2(F1_code, x, Q2);
                _F2 = SF_PDF->xfxQ2(F2_code, x, Q2);
                _F3 = SF_PDF->xfxQ2(F3_code, x, Q2);
            }

            break;
        }
        case 3: // CKMT
        {
            // Check if we're below Q0
            double Q2_eval = Q2;
            if (Q2 <= SQ(config.CKMT.Q0)) {
                Q2_eval = SQ(config.CKMT.Q0);
            }

            double f1 = SF_PDF->xfxQ2(F1_code, x, Q2_eval);
            double f2 = SF_PDF->xfxQ2(F2_code, x, Q2_eval);
            double f3 = SF_PDF->xfxQ2(F3_code, x, Q2_eval);

            if (Q2 <= SQ(config.CKMT.Q0)) {  // if below Q0, apply CKMT

                double chi = RescalingVariable(x, Q2);
                double chi_Q0 = RescalingVariable(x, SQ(config.CKMT.Q0));
                if ( (chi >= 1.0) || (config.sf_type == SFType::bottom) || (config.sf_type == SFType::top) ) {
                    _F1 = 0.0;
                    _F2 = 0.0;
                    _F3 = 0.0;
                } else {
                    double _F2_CKMT = F2_CKMT(chi, Q2);
                    double _F3_CKMT = F3_CKMT(chi, Q2);
                    double _F1_CKMT = F1_CKMT(_F2_CKMT, chi, Q2);

                    double _F2_CKMT_Q0 = F2_CKMT(chi_Q0, SQ(config.CKMT.Q0));
                    double _F3_CKMT_Q0 = F3_CKMT(chi_Q0, SQ(config.CKMT.Q0));
                    double _F1_CKMT_Q0 = F1_CKMT(_F2_CKMT_Q0, chi_Q0, SQ(config.CKMT.Q0));

                    _F1 = _F1_CKMT * (f1 / _F1_CKMT_Q0);
                    _F2 = _F2_CKMT * (f2 / _F2_CKMT_Q0);
                    _F3 = _F3_CKMT * (f3 / _F3_CKMT_Q0);
                }
            } else {
                _F1 = f1;
                _F2 = f2;
                _F3 = f3;
            }
            break;
        }
        case 4: // PCAC + CKMT
        {
            // Check if we're below Q0
            double Q2_eval = Q2;
            if (Q2 <= SQ(config.CKMT.Q0)) {
                Q2_eval = SQ(config.CKMT.Q0);
            }

            double f1 = SF_PDF->xfxQ2(F1_code, x, Q2_eval);
            double f2 = SF_PDF->xfxQ2(F2_code, x, Q2_eval);
            double f3 = SF_PDF->xfxQ2(F3_code, x, Q2_eval);

            if (Q2 <= SQ(config.CKMT.Q0)) {
                double chi = RescalingVariable(x, Q2);
                double chi_Q0 = RescalingVariable(x, SQ(config.CKMT.Q0));

                if ( (chi >= 1.0) || (config.sf_type == SFType::bottom) || (config.sf_type == SFType::top) ) {
                    _F1 = 0.0;
                    _F2 = 0.0;
                    _F3 = 0.0;
                } else {
                    double _F2_CKMT    = F2_CKMT(chi, Q2);
                    double _F3_CKMT    = F3_CKMT(chi, Q2);
                    double _F1_CKMT    = F1_CKMT(_F2_CKMT, chi, Q2);

                    double _F2_CKMT_Q0 = F2_CKMT(chi_Q0, SQ(config.CKMT.Q0));
                    double _F3_CKMT_Q0 = F3_CKMT(chi_Q0, SQ(config.CKMT.Q0));
                    double _F1_CKMT_Q0 = F1_CKMT(_F2_CKMT_Q0, chi_Q0, SQ(config.CKMT.Q0));

                    double _F2_PCAC = F2_PCAC(chi, Q2);
                    double _F2_PCAC_Q0 = F2_PCAC(chi_Q0, SQ(config.CKMT.Q0));

                    _F1 = _F1_CKMT * (f1 / _F1_CKMT_Q0);
                    _F2 = (_F2_CKMT + _F2_PCAC) * (f2 / (_F2_CKMT_Q0 + _F2_PCAC_Q0 ));
                    _F3 = _F3_CKMT * (f3 / _F3_CKMT_Q0);
                }
            } else {
                _F1 = f1;
                _F2 = f2;
                _F3 = f3;
            }
            break;
        }
    }

    return {_F1, _F2, _F3, _FL};
}

void StructureFunction::BuildGrids(string outpath) {
    const int Nx = config.SF.Nx;
    const int NQ2 = config.SF.NQ2;

    std::vector<double> x_arr;
    std::vector<double> Q2_arr;

    // Get the coefficients for parton calculation
    if (config.SF.mass_scheme == "parton") {
        GetCoefficients();
    }

    // Load the SF from LHAPDF
    if (mode == 2) {
        SF_PDF = config.Get_LHAPDF_SF(1); // Loads the plain APFEL structure functions
    } else if (mode == 3) {
        SF_PDF = config.Get_LHAPDF_SF(2); // Loads the TMC structure functions
    } else if (mode == 4) {
        SF_PDF = config.Get_LHAPDF_SF(2); // Loads the TMC structure functions
    }

    // Step sizes in log space, we do a slightly modified thing for the x spacing to get more points near x=1
    double d_log_Q2 = std::abs( std::log10(config.SF.Q2min) - std::log10(config.SF.Q2max) ) / (NQ2-1);
    double d_log_x  = std::abs( std::log10(config.SF.xmin) - std::log10(0.1)  ) / (Nx-60);

    // setup grid stuff
    f1_grid_fn = outpath + "/F1_"+config.current+"_"+config.projectile+"_"+config.target+"_"+config.sf_type_string+outsuffix;//+".grid";
    f2_grid_fn = outpath + "/F2_"+config.current+"_"+config.projectile+"_"+config.target+"_"+config.sf_type_string+outsuffix;//+".grid";
    f3_grid_fn = outpath + "/F3_"+config.current+"_"+config.projectile+"_"+config.target+"_"+config.sf_type_string+outsuffix;//+".grid";
    fL_grid_fn = outpath + "/FL_"+config.current+"_"+config.projectile+"_"+config.target+"_"+config.sf_type_string+outsuffix;//+".grid";

    std::ofstream f1_outfile;
    f1_outfile.open(f1_grid_fn + ".grid");
    std::ofstream f2_outfile;
    f2_outfile.open(f2_grid_fn + ".grid");
    std::ofstream f3_outfile;
    f3_outfile.open(f3_grid_fn + ".grid");
    std::ofstream fL_outfile;
    fL_outfile.open(fL_grid_fn + ".grid");

    // Get the Q2 and x values
    // size_t N_samples = Nx * NQ2;
    for (int i = 0; i < NQ2; i++) {
        double log_Q2 = std::log10(config.SF.Q2min) + i * d_log_Q2;
        Q2_arr.push_back(log_Q2);
        f1_outfile << log_Q2 << " "; f2_outfile << log_Q2 << " "; f3_outfile << log_Q2 << " "; fL_outfile << log_Q2 << " ";
    }
    f1_outfile << "\n"; f2_outfile << "\n"; f3_outfile << "\n"; fL_outfile << "\n";
    // for (int j = 0; j < Nx; j++) {
    //     double log_x = std::log10(config.SF.xmin) + j * d_log_x;
    //     x_arr.push_back(log_x);
    //     // std::cout << log_x << std::endl;
    //     f1_outfile << log_x << " "; f2_outfile << log_x << " "; f3_outfile << log_x << " "; fL_outfile << log_x << " ";
    // }

    for (int j = 0; j < Nx - 60; j++) {
        double log_x = std::log10(config.SF.xmin) + j * d_log_x;
        x_arr.push_back(log_x);
        // std::cout << log_x << std::endl;
        f1_outfile << log_x << " "; f2_outfile << log_x << " "; f3_outfile << log_x << " "; fL_outfile << log_x << " ";
    }
    // std::cout << std::endl;
    for (int j = 0; j < 40; j++) {
        double log_x =  std::log10(0.1 + j * (0.925 - 0.1)/(40));
        x_arr.push_back(log_x);
        // std::cout << log_x << std::endl;
        f1_outfile << log_x << " "; f2_outfile << log_x << " "; f3_outfile << log_x << " "; fL_outfile << log_x << " ";
    }
    // std::cout << std::endl;
    for (int j = 0; j < 20; j++) {
        double log_x =  std::log10(0.925 + j * (1.0 - 0.925)/(20-1));
        x_arr.push_back(log_x);
        // std::cout << log_x << std::endl;
        f1_outfile << log_x << " "; f2_outfile << log_x << " "; f3_outfile << log_x << " "; fL_outfile << log_x << " ";
    }
    f1_outfile << "\n"; f2_outfile << "\n"; f3_outfile << "\n"; fL_outfile << "\n";

    // Collect SF values and write grids
    for (int i = 0; i < NQ2; i++) {
        double log_Q2 = Q2_arr.at(i);
        double Q2 = std::pow(10.0, log_Q2);

        // if ((Q2 >= 1e2) && config.SF.dynamic_small_x && !config.SF.enable_small_x) {
        //     config.SF.enable_small_x = true;
        //     InitializeAPFEL();
        // } else if ((Q2 < 1e2) && config.SF.dynamic_small_x && config.SF.enable_small_x) {
        //     config.SF.enable_small_x = false;
        //     InitializeAPFEL();
        // }

        if (config.general.debug) {
            std::cout << "Q2 = " << Q2 << std::endl;
        }

        if ( (config.SF.mass_scheme != "parton") && (mode == 1)) {
            Set_Q_APFEL(std::sqrt(Q2));
        }

        double pf = 1.;  // prefactor for getting the NC right
        if ( (config.current == "NC") && (mode == 1)) {
            pf = 2./pow( Q2/(Q2 + pow(APFEL::GetZMass(),2))/4/APFEL::GetSin2ThetaW()/(1-APFEL::GetSin2ThetaW()),2);
        }

        for (int j = 0; j < Nx; j++) {
            double log_x = x_arr.at(j);
            double x = std::pow(10.0, log_x);

            double _F1, _F2, _F3, _FL;
            std::tie(_F1, _F2, _F3, _FL) = EvaluateSFs(x, Q2);

            if(!std::isfinite(_F1)) {
                // std::cerr << "F1 Infinite! Q2 = " << Q2 << ", x = " << x << ". Setting to zero." << std::endl;
                _F1 = 0.0;
            }
            if(!std::isfinite(_F2)) {
                // std::cerr << "F2 Infinite! Q2 = " << Q2 << ", x = " << x << ". Setting to zero." << std::endl;
                _F2 = 0.0;
            }
            if(!std::isfinite(_F3)) {
                // std::cerr << "F3 Infinite! Q2 = " << Q2 << ", x = " << x << ". Setting to zero." << std::endl;
                _F3 = 0.0;
            }
            if(!std::isfinite(_FL)) {
                _FL = 0.0;
            }
            // std::cout << x << "," << Q2 << ": " << _F2 << std::endl;

            // Write grid data
            f1_outfile << pf*_F1; f2_outfile << pf*_F2; f3_outfile << pf*_F3; fL_outfile << pf*_FL;
            if (j < Nx-1) {
              f1_outfile << ","; f2_outfile << ","; f3_outfile << ","; fL_outfile << ",";
            } else {
              f1_outfile << "\n"; f2_outfile << "\n"; f3_outfile << "\n"; fL_outfile << "\n";
            }
        }
    }
}

}