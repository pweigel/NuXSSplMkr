#include "cross_section.h"

namespace nuxssplmkr {

CrossSection::CrossSection(Configuration config) {
    sf_info = config.sf_info;
    pc = new nuxssplmkr::PhysConst();
    M_iso = 0.5*(pc->proton_mass + pc->neutron_mass);
}

void CrossSection::Load_Structure_Functions(string sf1_path, string sf2_path, string sf3_path) {
    Load_F1(sf1_path);
    Load_F2(sf2_path);
    Load_F3(sf3_path);
}

void CrossSection::Load_F1(string path) {
    if (F1_loaded) { // If we already loaded a spline, delete the old one
        F1 = photospline::splinetable<>();
    }
    F1.read_fits(path);
    std::cout << "Loaded F1!" << std::endl;
    // TODO: some sort of check that this worked
    F1_loaded = true;
}

void CrossSection::Load_F2(string path) {
    if (F2_loaded) { // If we already loaded a spline, delete the old one
        F2 = photospline::splinetable<>();
    }
    F2.read_fits(path);
    std::cout << "Loaded F2!" << std::endl;
    F2_loaded = true;
}

void CrossSection::Load_F3(string path) {
    if (F3_loaded) { // If we already loaded a spline, delete the old one
        F3 = photospline::splinetable<>();
    }
    F3.read_fits(path);
    std::cout << "Loaded F3!" << std::endl;
    F3_loaded = true;
}

void CrossSection::Set_Neutrino_Energy(double E) {
    ENU = E;
}

double CrossSection::_ds_dxdy(double* k) {
    double x = std::exp(k[0]);
    double y = std::exp(k[1]);
    // Jacobian = x * y, needed because we're integrating over log space
    return x * y * ds_dxdy(x, y);
}

double CrossSection::ds_dxdy(double E, double x, double y) {
    Set_Neutrino_Energy(E);
    return ds_dxdy(x, y);
}

double CrossSection::ds_dxdy(double x, double y) {
    double MW2 = sf_info.M_boson2 * SQ(pc->GeV); // TODO: This should happen where M_boson2 is?

    double s = 2 * M_iso * ENU;
    double Q2 = s * x * y;

    double prefactor = SQ(pc->GF) / (2 * M_PI * x); 
    double propagator = SQ( MW2 / (Q2 + MW2) );
    double jacobian = s * x; // from d2s/dxdQ2 --> d2s/dxdy

    // Constraints
    if (Q2 / SQ(pc->GeV) < 1) {
        return 0.0;
    } 
    else if (Q2 / SQ(pc->GeV) > sf_info.Q2max) {
        Q2 = sf_info.Q2max * SQ(pc->GeV); // freeze
    }
    
    std::array<double, 2> pt{{std::log10(Q2 / SQ(pc->GeV)), std::log10(x)}};

    std::array<int, 2> F1_splc;
    std::array<int, 2> F2_splc;
    std::array<int, 2> F3_splc;

    F1.searchcenters(pt.data(), F1_splc.data());
    F2.searchcenters(pt.data(), F2_splc.data());
    F3.searchcenters(pt.data(), F3_splc.data());

    double F1_val = F1.ndsplineeval(pt.data(), F1_splc.data(), 0);
    double F2_val = F2.ndsplineeval(pt.data(), F2_splc.data(), 0);
    double F3_val = F3.ndsplineeval(pt.data(), F3_splc.data(), 0);
    
    double term1 = y * ( y * x ) * F1_val;
    double term2 = ( 1 - y ) * F2_val;
    double term3 = ( x*y*(1-y/2) ) * F3_val;
    
    double xs = prefactor * jacobian * propagator * (term1 + term2 + term3);
    return xs / SQ(pc->cm); // TODO: Unit conversion outside of this function?
}

double CrossSection::_ds_dy(double k) {
    double x = std::pow(10, k);
    return x * ds_dxdy(x, _kernel_y);
}

double CrossSection::ds_dxdy_TMC() {
    return 0.;
}

double CrossSection::ds_dy_TMC() {
    return 0.;
}

double CrossSection::ds_dy(double E, double y) {
    Set_Neutrino_Energy(E);
    _kernel_y = y;

    gsl_integration_workspace * w = gsl_integration_workspace_alloc(5000);
    double result, error;

    gsl_function F;
    F.function = &KernelHelper<CrossSection, &CrossSection::_ds_dy>;
    F.params = this;
    
    gsl_integration_qag ( &F, log(1.e-8), log(1.), 0, 1.e-5, 5000, 6, w, &result, &error);
    gsl_integration_workspace_free(w);

    return result;
}

double CrossSection::TotalXS(double E){
    Set_Neutrino_Energy(E);

    double res,err;
    const unsigned long dim = 2; int calls = 50000;

    // integrating on the log of x and y
    double xl[dim] = { log(1.e-8), log(1.e-8) };
    double xu[dim] = { log(1.)   , log(1.)    };

    gsl_rng_env_setup ();
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *r = gsl_rng_alloc (T);

    gsl_monte_function F = { &KernelHelper<CrossSection, &CrossSection::_ds_dxdy>, dim, this};
    gsl_monte_vegas_state *s_vegas = gsl_monte_vegas_alloc (dim);

    gsl_monte_vegas_integrate (&F, xl, xu, dim, 10000, r, s_vegas, &res, &err);

    do
    {
        gsl_monte_vegas_integrate (&F, xl, xu, dim, calls, r, s_vegas, &res, &err);
    }
    while (fabs (gsl_monte_vegas_chisq (s_vegas) - 1.0) > 0.5 );

    gsl_monte_vegas_free (s_vegas);
    gsl_rng_free (r);

    return res;
}

}