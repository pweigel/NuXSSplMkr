#include "structure_function.h"

namespace nuxssplmkr {

StructureFunction::StructureFunction(Configuration &_config)
    : config(_config)
{
    pc = new PhysConst();

    GF2   = SQ(pc->GF);
    M_iso = 0.5*(pc->proton_mass + pc->neutron_mass);

    // Calculate fundamental constants
    // s_w = config.Sin2ThW;
    // Lu2 = ( 1. - (4./3.)*s_w) * ( 1. - (4./3.)*s_w);
    // Ld2 = (-1. + (2./3.)*s_w) * (-1. + (2./3.)*s_w);
    // Ru2 = (    - (4./3.)*s_w) * (    - (4./3.)*s_w);
    // Rd2 = (      (2./3.)*s_w) * (      (2./3.)*s_w);
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
    if (config.SF.disable_top == true) { // TODO: a better way of doing this. This should only be used for CSMS I think.
        std::cout << "WARNING: Top mass set to m_b + 0.1!" << std::endl;
        APFEL::SetPoleMasses(config.pdf.pdf_quark_masses[4], config.pdf.pdf_quark_masses[5], config.pdf.pdf_quark_masses[5]+0.1);
    } else {
        APFEL::SetPoleMasses(config.pdf.pdf_quark_masses[4], config.pdf.pdf_quark_masses[5], config.pdf.pdf_quark_masses[6]);
    }

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
    //APFEL::SetFFNS(3);
    //APFEL::SetTheory("QUniD");
    //APFEL::SetTheory("QED");
    //APFEL::SetTheory("QUniD");
    //APFEL::EnableLeptonEvolution(true);
    //APFEL::SetTauMass(1e10);
    //APFEL::SetPDFEvolution("exactalpha");
    APFEL::SetNumberOfGrids(3);
    APFEL::SetGridParameters(1, 90, 3, config.SF.xmin);
    APFEL::SetGridParameters(2, 80, 5, 1e-1);
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

double StructureFunction::NachtmannR(double x, double Q2){
    double m = M_iso / (pc->GeV);
    return sqrt(1. + 4.*x*x*SQ(m)/Q2);
}

double StructureFunction::NachtmannXi(double x, double Q2){
    double denominator = 1. + NachtmannR(x, Q2);
    return 2.*x / denominator;
}

template<class T,double (T::*f)(double,double), int n, int m>
double StructureFunction::HGeneric(double xi, double Q2){
    // Integrate T::f from xi to 1 at Q2
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
    double result, error;
    
    gsl_function F;
    Set_Q_APFEL(std::sqrt(Q2));
    _kernel_Q2 = Q2;
    F.function = &HK<T, f, n, m>;
    F.params = this;
    // double _integrate_xmax = 1.0;
    gsl_integration_qags ( &F, xi, 0.99, 0, 1.e-5, 5000, w, &result, &error);
    gsl_integration_workspace_free(w);

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
    double r = NachtmannR(x, Q2);
    double m = M_iso / (pc->GeV);
    double term1 = (x / (xi * r)) * F1(xi, Q2);
    double term2 = SQ(m * x / r) / Q2 * H2(xi, Q2);
    double term3 = 2.0 * SQ(SQ(m) * x / (Q2 * r)) * (x / r) * G2(xi, Q2);
    return (term1 + term2 + term3);
}

double StructureFunction::F2_TMC(double x, double Q2) {
    double xi = NachtmannXi(x, Q2);
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
    double r = NachtmannR(x, Q2);
    double m = M_iso / (pc->GeV);

    // double term1 = (x / (xi * SQ(r))) * F3(xi, Q2);
    // double term2 = 2.0 * SQ(m * x / r) / (Q2 * r) * H3(xi, Q2)
    double g2 = G2(xi, Q2);
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
    // reno paper
    double delta = CKMT_Delta(Q2);
    double n = CKMT_n(Q2);
    // double term1 = A * pow(x, -1.0 * delta) *  pow((1.0 - x), n+4);
    // double term2 = pow(Q2 / (Q2 + config.CKMT.a), 1 + delta);
    // double term3 = B * pow(x, 1.0 - config.CKMT.AlphaR) * pow(1.0 - x, n) * pow(Q2 / (Q2 + config.CKMT.b), config.CKMT.AlphaR);
    // double term4 = (1.0 + f * (1.0 - x));
    // return (term1 * term2 + term3 * term4);
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

// double StructureFunction::FL_PCAC(double x, double Q2) {
//     return 0.0;
// }

// double StructureFunction::F2_PCAC(double x, double Q2) {
//     return 0.0;
// }

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
    for ( double log_Q2 = std::log10(config.SF.Q2min) - d_log_Q2_knot; log_Q2 <= std::log10(config.SF.Q2max) + d_log_Q2_knot; log_Q2 += d_log_Q2_knot ) {
        double knot = log_Q2;
        Q2_knots.push_back(knot);
    }

    for ( double log_x = std::log10(config.SF.xmin) - d_log_x_knot; log_x <= 1 + d_log_x_knot; log_x += d_log_x_knot ) {
        double knot = log_x;
        x_knots.push_back(knot);
    }

    // Collect SF values
    std::deque<std::pair<double,std::array<unsigned int, 2>>> F1_spline_data;
    std::deque<std::pair<double,std::array<unsigned int, 2>>> F2_spline_data;
    std::deque<std::pair<double,std::array<unsigned int, 2>>> F3_spline_data;
    for (unsigned int Q2i = 0; Q2i < NQ2; Q2i++) {
        double log_Q2 = std::log10(config.SF.Q2min) + Q2i * d_log_Q2;
        double Q2 = std::pow(10.0, log_Q2);

        if (config.SF.mass_scheme != "parton") {
            Set_Q_APFEL(std::sqrt(Q2));
        }

        for (unsigned int xi = 0; xi < Nx; xi++) {
            double log_x = std::log10(config.SF.xmin) + xi * d_log_x;
            double x = std::pow(10.0, log_x);

            // Do checks here
           // if ( Q2*(1/z-1)+mass_nucl*mass_nucl <= TMath::Power(mass_nucl+mPDFQrk[TMath::Abs(pdg_fq)],2) ) { sf_stream << 0. << "  "; continue; }
            // if (Q2 * (1/x - 1) + 0.93 * 0.93 <= )

            double _FL, _F1, _F2, _F3;
            _FL = FL(x, Q2);
            _F2 = F2(x, Q2);
            _F1 = (_F2 - _FL) / (2. * x);
            _F3 = F3(x, Q2);
            if (config.SF.enable_CKMT) {
                double CKMT_Q20 = SQ(config.CKMT.Q0);
                if (Q2 < CKMT_Q20) {
                    double _F2_CKMT    = F2_CKMT(x, Q2);
                    double _F3_CKMT    = F3_CKMT(x, Q2);
                    double _F2_Q0      = F2(x, CKMT_Q20);
                    double _F3_Q0      = F3(x, CKMT_Q20);
                    double _F2_CKMT_Q0 = F2_CKMT(x, CKMT_Q20);
                    double _F3_CKMT_Q0 = F3_CKMT(x, CKMT_Q20);

                    // if (config.SF.enable_PCAC) {
                        
                    // }

                    // R parameterization from Whitlow et al: Phys. Lett. B 250, 193 (1990)
                    // R is used to calculate F1 from F2
                    double big_theta = 1.0 + 12.0 * (Q2 / (Q2 + 1.0)) * (SQ(0.125) / (SQ(0.125) + SQ(x)));
                    // double b1 = 0.635;
                    // double b2 = 0.5747;
                    // double b3 = -0.3534;
                    double _R = 0.635 / log(Q2/SQ(0.2)) * big_theta + 0.5747 / Q2 - 0.3534 / (Q2 + SQ(0.3));
                    
                    _F2 = _F2_CKMT * (_F2_Q0 / _F2_CKMT_Q0);
                    _F1 = _F2 * (1.0 + 4.0 * SQ(M_iso * x)) / (2.0 * x * (_R + 1.0));
                    _F3 = _F3_CKMT * (_F3_Q0 / _F3_CKMT_Q0);
                }
            }

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

            F1_spline_data.push_back(std::make_pair(_F1, std::array<unsigned int, 2>{Q2i, xi}));
            F2_spline_data.push_back(std::make_pair(_F2, std::array<unsigned int, 2>{Q2i, xi}));
            F3_spline_data.push_back(std::make_pair(_F3, std::array<unsigned int, 2>{Q2i, xi}));
        }
    }

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

    double smooth = 1e-10;

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

void StructureFunction::BuildGrids(string outpath) {
    const unsigned int Nx = config.SF.Nx;
    const unsigned int NQ2 = config.SF.NQ2;

    std::vector<double> x_arr;
    std::vector<double> Q2_arr;

    // Get the coefficients for parton calculation
    GetCoefficients();

    std::ofstream F1_file;
    std::ofstream F2_file;
    std::ofstream F3_file;
    F1_file.open(outpath + "/F1_"+config.projectile+"_"+config.target+"_"+config.sf_type_string+".grid");
    F2_file.open(outpath + "/F2_"+config.projectile+"_"+config.target+"_"+config.sf_type_string+".grid");
    F3_file.open(outpath + "/F3_"+config.projectile+"_"+config.target+"_"+config.sf_type_string+".grid");

    // Step sizes in log space
    double d_log_Q2 = std::abs( std::log10(config.SF.Q2min) - std::log10(config.SF.Q2max) ) / NQ2;
    double d_log_x  = std::abs( std::log10(config.SF.xmin)  - std::log10(config.SF.xmax)  ) / Nx;

    for (unsigned int Q2i = 0; Q2i < NQ2; Q2i++) {
        double log_Q2 = std::log10(config.SF.Q2min) + (0.5 + Q2i) * d_log_Q2;
        double Q2 = std::pow(10.0, log_Q2);

        Set_Q_APFEL(std::sqrt(Q2));

        for (unsigned int xi = 0; xi < Nx; xi++) {
            double log_x = std::log10(config.SF.xmin) + (0.5 + xi) * d_log_x;
            double x = std::pow(10.0, log_x);

            // Do checks here
           // if ( Q2*(1/z-1)+mass_nucl*mass_nucl <= TMath::Power(mass_nucl+mPDFQrk[TMath::Abs(pdg_fq)],2) ) { sf_stream << 0. << "  "; continue; }
            // if (Q2 * (1/x - 1) + 0.93 * 0.93 <= )
            double _FL = FL(x, Q2); 
            double _F2 = F2(x, Q2);
            // calculate F1 from FL, F2 instead of calling F1(x, Q2), which recomputes
            double _F1 = (_F2 - _FL) / (2. * x);
            double _F3 = F3(x, Q2);

            if(!std::isfinite(_F1)) {
                std::cerr << "F1 Infinite! Q2 = " << Q2 << ", x = " << x << ". Setting to zero." << std::endl;
                _F1 = 0.0;
            } else if (_F1 < 0) {
                _F1 = 0.0;
            }
            if(!std::isfinite(_F2)) {
                std::cerr << "F2 Infinite! Q2 = " << Q2 << ", x = " << x << ". Setting to zero." << std::endl;
                _F2 = 0.0;
            } else if (_F2 < 0) {
                _F2 = 0.0;
            }
            if(!std::isfinite(_F3)) {
                std::cerr << "F3 Infinite! Q2 = " << Q2 << ", x = " << x << ". Setting to zero." << std::endl;
                _F3 = 0.0;
            } else if (_F3 < 0) {
                _F3 = 0.0;
            }

            F1_file << log_Q2 << "," << log_x << "," << _F1 << "\n";
            F2_file << log_Q2 << "," << log_x << "," << _F2 << "\n";
            F3_file << log_Q2 << "," << log_x << "," << _F3 << "\n";
        }
    }
}

// double StructureFunction::Evaluate(double Q2, double x, double y){
//     // only evaluates central values

//     LHAPDF::GridPDF* grid_central = dynamic_cast<LHAPDF::GridPDF*>(config.pdf.pdf);
//     string xt = "nearest";
//     grid_central -> setExtrapolator(xt);

//     map<int,double> xq_arr;
//     for ( int p : partons ){
//       xq_arr[p] = grid_central -> xfxQ2(p, x, Q2 / (pc->GeV2));
//     }
    
//     return SigR_Nu_LO(x, y, xq_arr);
// }

// double StructureFunction::SigR_Nu_LO(double x, double y, map<int,double> xq_arr){
//     double k = 0.;
//     d_lepton = SQ(M_lepton)/(2.*M_iso*ENU);
//     double y_p = (1. - d_lepton / x) + (1.- d_lepton/x - y) * (1. - y);
//     double y_m = (1. - d_lepton / x) - (1.- d_lepton/x - y) * (1. - y);
//     double a = y_p + CP_factor*y_m;
//     double b = y_p - CP_factor*y_m;

//       map<int,double> SigRcoef;

//       // Coefficients for CC
//     SigRcoef[1]  =    a ;
//     SigRcoef[-1] =    b ;
//     SigRcoef[2]  =    a ;
//     SigRcoef[-2] =    b ;
//     SigRcoef[3]  = 2.*a ;
//     SigRcoef[-3] =   0. ;
//     SigRcoef[4]  =   0. ;
//     SigRcoef[-4] = 2.*b ;
//     SigRcoef[5]  = 2.*a ;
//     SigRcoef[-5] =   0. ;
//     SigRcoef[21] =   0. ;

//     if (CP_factor < 0 ){
//         //fixes for antineutrinos
//         SigRcoef[3]   = 0. ;
//         SigRcoef[-3]  = 2.*b ;
//         SigRcoef[4]   = 2.*a ;
//         SigRcoef[-4]  = 0. ;
//         SigRcoef[5]   = 0. ;
//         SigRcoef[-5]  = 2.*b ;
//     }

//     for( int p : partons ) {
//         k += SigRcoef[p]*xq_arr[p];
//     }

//     return k;
// }

}