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

    CP_factor = CPFactorMap.at(sf_info.neutrino_type);

    if (sf_info.current == CC) {
        M_boson2 = SQ(pc->Wboson_mass);
    } else if (sf_info.current == NC) {
        M_boson2 = SQ(pc->Zboson_mass);
    } else {
        throw std::runtime_error("Unidentified current specified!");
    }

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
    APFEL::SetGridParameters(1, 90, 3, sf_info.xmin);
    APFEL::SetGridParameters(2, 50, 5, 2e-1);
    APFEL::SetGridParameters(3, 40, 5, 8e-1);
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

void StructureFunction::GetCoefficients() {
    int down_type[3] = {1, 3, 5};  // d s b
    int up_type[3] = {2, 4, 6};  // u c t
    // TODO: Maybe do isospin symmetry here?

    if (sf_info.neutrino_type == neutrino) {
        for (int i : down_type) {
            F2coef[i] = 1.;
            F3coef[i] = 1.;
        }
        for (int i: up_type) {
            F2coef[-i] = 1.;
            F2coef[-i] = -1.;
        }
    } else if (sf_info.neutrino_type == antineutrino) {
        for (int i : up_type) {
            F2coef[i] = 1.;
            F3coef[i] = 1.;
        }
        for (int i: down_type) {
            F2coef[-i] = 1.;
            F2coef[-i] = -1.;
        }
    } else {
        throw std::runtime_error("Unidentified neutrino type!");
    }

    // gluons
    F2coef[21] = 0.;
    F3coef[21] = 0.; 
}

double StructureFunction::F1(double x, double Q2) {
    // LO for now
    if ( (sf_info.perturbative_order == LO) && (!sf_info.Use_APFEL_LO) ) {
        return F2(x, Q2) / (2. * x);
    } else {
        return (F2(x, Q2) - FL(x, Q2)) / (2. * x);
    }
}

double StructureFunction::F2(double x, double Q2) {
    if ( (sf_info.perturbative_order == LO) && (!sf_info.Use_APFEL_LO) ) {
        auto s = PDFExtract(x, Q2);
        
        // TODO: Figure out a better way to do this (isospin symmetry for p<->n)
        if(sf_info.target_type == neutron) {
            auto q1 = s[1];
            auto q2 = s[2];
            s[1] = q2;
            s[2] = q1;
        }
        return F2_LO(s);
    } else {
        return APFEL::F2total(x);
    }
}

double StructureFunction::FL(double x, double Q2) {
    if ( (sf_info.perturbative_order == LO) && (!sf_info.Use_APFEL_LO) ) {
        return F1(x, Q2) - 2 * x * F2(x, Q2);
    } else {
        return APFEL::FLtotal(x);
    }
}

double StructureFunction::F2_LO(map<int, double>& xq_arr) {
    /* 
    Example: Neutrino
    F2 = 2x(d + s + b + ubar + cbar + tbar)
    */
    double k = 0.;
    for( int p : partons ) {
        k += F2coef[p] * xq_arr[p];
    }	
    return 2 * k;
}

double StructureFunction::xF3(double x, double Q2) {
    // only true at LO

    if ( (sf_info.perturbative_order == LO) && (!sf_info.Use_APFEL_LO) ) {
        auto s = PDFExtract(x, Q2);

        // TODO: Figure out a better way to do this (isospin symmetry for p<->n)
        if(sf_info.target_type == neutron) {
            auto q1 = s[1];
            auto q2 = s[2];
            s[1] = q2;
            s[2] = q1;
        }
        return xF3_LO(s);
    } else {
        return APFEL::F3total(x);
    }
}

double StructureFunction::xF3_LO(map<int, double>& xq_arr) {
    /* Example: Neutrino
    xF3 = 2x(d + s + b - ubar - cbar - tbar)
    */
    double k=0.;
    for( int p : partons ) {        
        k += F3coef[p]*xq_arr[p];
    }	
    return 2 * k;
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
      xq_arr[p] = grid_central -> xfxQ2(p, x, Q2/(pc->GeV2));
    }
    return xq_arr;
}

void StructureFunction::BuildGrids(string outpath) {
    const unsigned int Nx = sf_info.Nx;
    const unsigned int NQ2 = sf_info.NQ2;

    std::vector<double> x_arr;
    std::vector<double> Q2_arr;

    // Step sizes in log space
    double d_log_Q2 = std::abs( std::log10(sf_info.Q2min) - std::log10(sf_info.Q2max) ) / NQ2;
    double d_log_x  = std::abs( std::log10(sf_info.xmin)  - std::log10(sf_info.xmax)  ) / Nx;

    // Spline parameters
    const uint32_t dim = 2;
    std::vector<uint32_t> orders(dim, 2);

    unsigned int Nknots_Q2 = 100;
    unsigned int Nknots_x = 100;

    std::vector<double> Q2_knots;
    std::vector<double> x_knots;

    // Step sizes for knots in log space
    const double d_log_Q2_knot= std::abs(std::log10(sf_info.Q2max)- std::log10(sf_info.Q2min)) / (Nknots_Q2 - 1);
    const double d_log_x_knot = std::abs(std::log10(sf_info.xmax) - std::log10(sf_info.xmin) ) / (Nknots_x - 1);

    std::cout << "log_Q2min = " << std::log10(sf_info.Q2min) << ", log_Q2max = " << std::log10(sf_info.Q2max) << std::endl;
    std::cout << "log_xmin = " << std::log10(sf_info.xmin) << ", log_xmax = " << std::log10(sf_info.xmax) << std::endl;
    std::cout << "d_log_Q2 = " << d_log_Q2 << ", d_log_x = " << d_log_x << std::endl;

    // Get the Q2 and x values
    size_t N_samples = Nx * NQ2;
    for ( double log_Q2 = std::log10(sf_info.Q2min); log_Q2<std::log10(sf_info.Q2max); log_Q2 += d_log_Q2 ) {
        double powQ2 = std::pow( 10, log_Q2 + 0.5 * d_log_Q2 );
        double Q2 = log_Q2 + 0.5 * d_log_Q2;
        if (powQ2 > sf_info.Q2max) continue;
        Q2_arr.push_back(Q2);
    }

    for ( double log_x = std::log10(sf_info.xmin); log_x<std::log10(sf_info.xmax); log_x += d_log_x ) {
        double powx = std::pow( 10, log_x + 0.5*d_log_x );
        double x = log_x + 0.5*d_log_x;
        if ( powx > sf_info.xmax ) continue;
        x_arr.push_back(x);
    }

    // Get the knots for Q2 and x
    for ( double log_Q2 = std::log10(sf_info.Q2min) - 2*d_log_Q2_knot; log_Q2 <= std::log10(sf_info.Q2max); log_Q2 += d_log_Q2_knot ) {
        double knot = log_Q2 + 0.5 * d_log_Q2_knot;
        Q2_knots.push_back(knot);
    }

    for ( double log_x = std::log10(sf_info.xmin) - 2*d_log_x_knot; log_x <= 1; log_x += d_log_x_knot ) {
        double knot = log_x + 0.5 * d_log_x_knot;
        x_knots.push_back(knot);
    }

    // Collect SF values
    std::deque<std::pair<double,std::array<unsigned int, 2>>> F1_spline_data;
    std::deque<std::pair<double,std::array<unsigned int, 2>>> F2_spline_data;
    std::deque<std::pair<double,std::array<unsigned int, 2>>> F3_spline_data;
    for (unsigned int Q2i = 0; Q2i < NQ2; Q2i++) {
        double log_Q2 = std::log10(sf_info.Q2min) + (0.5 + Q2i) * d_log_Q2;
        double Q2 = std::pow(10.0, log_Q2);

        Set_Q_APFEL(std::sqrt(Q2));

        for (unsigned int xi = 0; xi < Nx; xi++) {
            double log_x = std::log10(sf_info.xmin) + (0.5 + xi) * d_log_x;
            double x = std::pow(10.0, log_x);

            double _FL = FL(x, Q2);
            double _F2 = F2(x, Q2);
            // calculate F1 from FL, F2 instead of calling F1(x, Q2), which recomputes
            double _F1 = (_F2 - _FL) / (2. * x);
            double _F3 = F3(x, Q2);

            if(!std::isfinite(_F1)) {
                std::cerr << "F1 Infinite! Q2 = " << Q2 << ", x = " << x << ". Setting to zero." << std::endl;
                _F1 = 0.0;
            }
            if(!std::isfinite(_F2)) {
                std::cerr << "F2 Infinite! Q2 = " << Q2 << ", x = " << x << ". Setting to zero." << std::endl;
                _F2 = 0.0;
            }
            if(!std::isfinite(_F3)) {
                std::cerr << "F2 Infinite! Q2 = " << Q2 << ", x = " << x << ". Setting to zero." << std::endl;
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
    F1_spline.fit(F1_data, F1_weights, std::vector<std::vector<double>>{Q2_arr, x_arr}, std::vector<uint32_t>{2,2}, 
                  {Q2_knots, x_knots}, {smooth, smooth}, {2, 2});
    if(std::isnan(*F1_spline.get_coefficients())){
				std::cerr << "F1 spline fit has failed!" << std::endl;
    }

    photospline::splinetable<> F2_spline;
    F2_spline.fit(F2_data, F2_weights, std::vector<std::vector<double>>{Q2_arr, x_arr}, std::vector<uint32_t>{2,2}, 
                  {Q2_knots, x_knots}, {smooth, smooth}, {2, 2});
    if(std::isnan(*F2_spline.get_coefficients())){
				std::cerr << "F2 spline fit has failed!" << std::endl;
    }

    photospline::splinetable<> F3_spline;
    F3_spline.fit(F3_data, F3_weights, std::vector<std::vector<double>>{Q2_arr, x_arr}, std::vector<uint32_t>{2,2}, 
                  {Q2_knots, x_knots}, {smooth, smooth}, {2, 2});
    if(std::isnan(*F3_spline.get_coefficients())){
				std::cerr << "F3 spline fit has failed!" << std::endl;
    }

    // Write splines
    F1_spline.write_fits(outpath + "/F1.fits");
    F2_spline.write_fits(outpath + "/F2.fits");
    F3_spline.write_fits(outpath + "/F3.fits");
}

double StructureFunction::ds_dxdy(double x, double y, double Q2) {
    double _FL = FL(x, Q2);
    double _F2 = F2(x, Q2);
    double _F3 = xF3(x, Q2) / x;
    // calculate manually so we don't recomopute w/ APFEL
    double _F1 = (_F2 - _FL) / (2. * x);

    // TODO: is this defined somewhere else
    double norm;
    if (sf_info.current == CC) {
        norm = 1;
    } else if (sf_info.current == NC) {
        norm = 1;  // TODO: Fix norm for apfel NC!
    } else {

    }

    double cp_factor;  // TODO: Is this defined somewhere else?
    if (sf_info.neutrino_type == neutrino) {
        cp_factor = 1;
    } else if (sf_info.neutrino_type == antineutrino) {
        cp_factor = -1;
    } else {

    }

    return x*y*y * _F1 + (1 - y) * _F2 + cp_factor * x * y * (1 - y/2) * _F3;
}

double StructureFunction::Evaluate(double Q2, double x, double y){
    // only evaluates central values

    LHAPDF::GridPDF* grid_central = dynamic_cast<LHAPDF::GridPDF*>(sf_info.pdf);
    string xt = "nearest";
    grid_central -> setExtrapolator(xt);

    map<int,double> xq_arr;
    for ( int p : partons ){
      xq_arr[p] = grid_central -> xfxQ2(p, x, Q2 / (pc->GeV2));
    }
    
    return SigR_Nu_LO(x, y, xq_arr);
}

double StructureFunction::SigR_Nu_LO(double x, double y, map<int,double> xq_arr){
	double k = 0.;
    d_lepton = SQ(M_lepton)/(2.*M_iso*ENU);
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

double StructureFunction::KernelXS(double * k){
    double x = exp(k[0]);
    double y = exp(k[1]);
    double s = 2.*M_iso*ENU + SQ(M_iso);
    double Q2 = ( s - SQ(M_iso) )*x*y;

    // std::cout << "Q2 = " << Q2 << std::endl;
    // std::cout << "Q2/SQ(pc->GeV) = " << Q2/SQ(pc->GeV) << std::endl;
    // bool cond = Q2/SQ(pc->GeV) < 0.6;
    // std::cout << "Q2/SQ(pc->GeV) < 5.0 = " << cond << std::endl;

    // if(Q2/SQ(pc->GeV) < 0.6){
    //     return 1.e-99;
    // }

    double denum    = SQ(1. + Q2/M_boson2);
    double norm     = GF2*M_iso*ENU/(2.*M_PI*denum);

    d_lepton = SQ(M_lepton)/(2.*M_iso*ENU);
    d_nucleon = M_iso / (2. * ENU);
    //   if(INT_TYPE==CC || IS_HNL==true){  // only CC, but if it's HNL then also NC
    //     //Following HEP PH 0407371 Eq. (7)
    //     double h = x*y + d_lepton;
    //     if((1. + x* d_nucleon) * h*h - (x+ d_lepton)*h + x * d_lepton > 0.){
    //         return 0.;
    //     }
    //   }
    double h = x*y + d_lepton;
    if((1. + x* d_nucleon) * h*h - (x+ d_lepton)*h + x * d_lepton > 0.){
        return 0.;
    }

    // x*y is the jacobian
    // std::cout << "x*y*norm*Evaluate(Q2, x, y) = " << x*y*norm*Evaluate(Q2, x, y) << std::endl;
    // std::cout << "Using function: LHAXS::KernelXS(double * k)" << std::endl;
    return x*y*norm*Evaluate(Q2, x, y);
}

double StructureFunction::TotalXS(){
    // TODO:
    // We used to do 10k warmup calls, 50k final calls
    // for the integrator. I changed this because it took forever (another issue).

    double res,err;
    const unsigned long dim = 2; int calls = 50000;

    // integrating on the log of x and y
    double xl[dim] = { log(1.e-7) , log(1.e-7) };
    double xu[dim] = { log(1.) , log(1.)};

    gsl_rng_env_setup ();
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *r = gsl_rng_alloc (T);

    gsl_monte_function F = { &KernelHelper<StructureFunction, &StructureFunction::KernelXS>, dim, this};
    gsl_monte_vegas_state *s_vegas = gsl_monte_vegas_alloc (dim);
    // std::cout << "Starting first integration... ";
    // TODO: Integration tests! -PW
    gsl_monte_vegas_integrate (&F, xl, xu, dim, 1000, r, s_vegas, &res, &err);
    // gsl_monte_vegas_integrate (&F, xl, xu, dim, 10000, r, s_vegas, &res, &err);
    // std::cout << " Done!" << std::endl;
    // std::cout << "Starting second integration... ";
    // do {
    //     gsl_monte_vegas_integrate (&F, xl, xu, dim, calls, r, s_vegas, &res, &err);
    // }
    // while (fabs (gsl_monte_vegas_chisq (s_vegas) - 1.0) > 0.5 );
    // std::cout << "Done!" << std::endl;

    gsl_monte_vegas_free (s_vegas);
    gsl_rng_free (r);

    return res;
}

void StructureFunction::Set_Lepton_Mass(double m) {
    M_lepton = m;
}

void StructureFunction::Set_Neutrino_Energy(double E) {
    ENU = E;
}

void StructureFunction::Set_Use_APFEL_LO(bool value) {
    sf_info.Use_APFEL_LO = value;
}

void StructureFunction::Set_Q_APFEL(double Q) {
    APFEL::SetAlphaQCDRef(sf_info.pdf->alphasQ(Q), Q);
    APFEL::ComputeStructureFunctionsAPFEL(Q, Q);
}

}