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
    F1.read_fits(path);
    std::cout << "Loaded F1!" << std::endl;
    // TODO: some sort of check that this worked
    F1_loaded = true;
}

void CrossSection::Load_F2(string path) {
    F2.read_fits(path);
    std::cout << "Loaded F2!" << std::endl;
    F2_loaded = true;
}

void CrossSection::Load_F3(string path) {
    F3.read_fits(path);
    std::cout << "Loaded F3!" << std::endl;
    F3_loaded = true;
}

void CrossSection::Set_Energy(double E) {
    ENU = E;
}

double CrossSection::ds_dxdy(double* k) {
    double x = std::pow(10, k[0]);
    double y = std::pow(10, k[1]);
  
    double s = (2.0 * M_iso * ENU + SQ(M_iso));
    double Q2 = ( s - SQ(M_iso) )*x*y;
    double prefactor_denom = SQ(1 + Q2 / (sf_info.M_boson2 * SQ(pc->GeV) ));
    double prefactor = SQ(pc->GF) * M_iso * ENU /(2.0 * M_PI * prefactor_denom);

    if (Q2 / SQ(pc->GeV) < 1e-2) {
        return 0.0;
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

    double term1 = y * ( x*y ) * F1_val;
    double term2 = ( 1 - y ) * F2_val;
    double term3 = ( x*y*(1-y/2) ) * F3_val;

    // x*y = Jacobian?
    return prefactor * (term1 + term2 + term3) / SQ(pc->cm);
    // return x * y * prefactor * (term1 + term2 + term3) / SQ(pc->cm);
    // return prefactor * ( term1*F1(x, Q2) + term2*F2(x, Q2) + term3*F3(x, Q2));
}

double CrossSection::ds_dxdy(double x, double y, double E) {
    double s = (2.0 * M_iso * E + SQ(M_iso));
    double Q2 = ( s - SQ(M_iso) )*x*y;
    double prefactor_denom = SQ(1 + Q2 / (sf_info.M_boson2 * SQ(pc->GeV) ));
    double prefactor = SQ(pc->GF) * M_iso * E /(2.0 * M_PI * prefactor_denom);

    if (Q2 / SQ(pc->GeV) < 1e-2) {
        return 0.0;
    }
    
    // std::cout << prefactor << " " << prefactor_denom << std::endl;

    std::array<double, 2> pt{{std::log10(Q2 / SQ(pc->GeV)), std::log10(x)}};
    // std::cout << std::log10(Q2 / SQ(pc->GeV)) << " " << std::log10(x) << std::endl;
	  std::array<int, 2> F1_splc;
	  std::array<int, 2> F2_splc;
	  std::array<int, 2> F3_splc;
    F1.searchcenters(pt.data(), F1_splc.data());
    F2.searchcenters(pt.data(), F2_splc.data());
    F3.searchcenters(pt.data(), F3_splc.data());

    double F1_val = F1.ndsplineeval(pt.data(), F1_splc.data(), 0);
    double F2_val = F2.ndsplineeval(pt.data(), F2_splc.data(), 0);
    double F3_val = F3.ndsplineeval(pt.data(), F3_splc.data(), 0);

    double term1 = y * ( x*y ) * F1_val;
    double term2 = ( 1 - y ) * F2_val;
    double term3 = ( x*y*(1-y/2) ) * F3_val;
    // x*y = Jacobian?
    return prefactor * (term1 + term2 + term3) / SQ(pc->cm);
    // return x * y * prefactor * (term1 + term2 + term3) / SQ(pc->cm);
    // return prefactor * ( term1*F1(x, Q2) + term2*F2(x, Q2) + term3*F3(x, Q2));
}

double CrossSection::ds_dy() {
    return 0.;
}

double CrossSection::ds_dxdy_TMC() {
    return 0.;
}

double CrossSection::ds_dy_TMC() {
    return 0.;
}

// double CrossSection::TrapezoidalTotalXS(double E) {
    
// }

double CrossSection::TotalXS(double E){
    // TODO:
    // We used to do 10k warmup calls, 50k final calls
    // for the integrator. I changed this because it took forever (another issue).
    Set_Energy(E);

    double res,err;
    const unsigned long dim = 2; int calls = 5000;

    // integrating on the log of x and y
    double xl[dim] = { std::log10(1.e-4), std::log10(1.e-4) };
    double xu[dim] = { std::log10(1.)   , std::log10(1.)    };

    gsl_rng_env_setup ();
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *r = gsl_rng_alloc (T);

    gsl_monte_function F = { &KernelHelper<CrossSection, &CrossSection::ds_dxdy>, dim, this};
    gsl_monte_vegas_state *s_vegas = gsl_monte_vegas_alloc (dim);
    gsl_monte_vegas_integrate (&F, xl, xu, dim, calls, r, s_vegas, &res, &err);
    gsl_monte_vegas_free (s_vegas);
    gsl_rng_free (r);

    return res;
}

}