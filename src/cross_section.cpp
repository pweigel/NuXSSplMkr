#include "NuXSSplMkr/cross_section.h"

namespace nuxssplmkr {

CrossSection::CrossSection(Configuration& _config, PhaseSpace& _ps)  
    : config(_config), ps(_ps)
{
    pc = new nuxssplmkr::PhysConst();
    rc_prefactor = (pc->alpha / (2.0 * M_PI));
    Set_Mode(config.mode);

    F1_code = std::stoi(config.Get_SF_Code("F1"));
    F2_code = std::stoi(config.Get_SF_Code("F2"));
    F3_code = std::stoi(config.Get_SF_Code("F3"));

    std::cout << "SF codes: " << std::endl;
    std::cout << F1_code << std::endl;
    std::cout << F2_code << std::endl;
    std::cout << F3_code << std::endl;

    SF_PDF = config.Get_LHAPDF_SF(mode);

    // TODO: not used yet, but needed for NC
    // double coef = pc->GF * pow(top_mass * pc->GeV, 2) / 8 / sqrt(2) / pow(M_PI, 2);
    // double rho  = 1 + 3 * coef * ( 1 + coef * ( 19 - 2 * pow(M_PI, 2) ) );
}

void CrossSection::Set_Mode(int _mode) {
    mode = _mode;

    switch (mode) {
        case 0: {insuffix = ""; outsuffix = ""; break;}
        case 1: {insuffix = ""; outsuffix = ""; break;}
        case 2: {insuffix = "_TMC"; outsuffix = "_TMC"; break;}
        case 3: {insuffix = "_CKMT"; outsuffix = "_CKMT"; break;}
        case 4: {insuffix = "_PCAC"; outsuffix = "_PCAC"; break;}
    }
}

void CrossSection::Set_Neutrino_Energy(double E) {
    ENU = E;
    integral_max_Q2 = min(config.XS.Q2max*SQ(pc->GeV), 2.0 * config.target_mass * ENU);
}

double CrossSection::ds_dxdy_kernel(double* k) {
    double x = std::exp(k[0]);
    double y = std::exp(k[1]);

    bool ps_valid = ps.Validate(ENU, x, y);
    if (!ps_valid) {
        return 0;
    }

    double ddxs = ds_dxdy(x, y);
    if (config.XS.enable_radiative_corrections) {
        double radiative_correction = rc_dsdxdy(ENU, x, y, ddxs);
        ddxs += radiative_correction;
        // double radiative_correction = rc_bardin(ENU, x, y);
        // ddxs *= (1.0 + radiative_correction);
    }

    // Jacobian = x * y, needed because we're integrating over log space
    return x * y * ddxs;
}

double CrossSection::ds_dxdy(double E, double x, double y) {
    Set_Neutrino_Energy(E);
    return ds_dxdy(x, y);
}


double CrossSection::ds_dxdQ2_kernel(double* k) {
    double x = std::exp(k[0]);
    double Q2 = std::exp(k[1]);

    bool ps_valid = ps.Validate_xQ2(ENU, x, Q2);
    if (!ps_valid) {
        return 0;
    }

    // Jacobian = x * Q2, needed because we're integrating over log space
    return x * Q2 * ds_dxdQ2(x, Q2);
}

double CrossSection::ds_dxdQ2(double E, double x, double Q2) {
    Set_Neutrino_Energy(E);
    return ds_dxdQ2(x, Q2);
}

double CrossSection::ds_dxdQ2(double x, double Q2) {
    double MW2 = config.constants.Mboson2 * SQ(pc->GeV); // TODO: This should happen where M_boson2 is?
    double M_target = config.target_mass;
    double M_l = config.lepton_mass;

    double s = 2.0 * M_target * ENU + SQ(M_target);
    double y = Q2 / ( (s - SQ(M_target)) * x);

    double prefactor = SQ(pc->GF) / (2 * M_PI * x); 
    double propagator = SQ( MW2 / (Q2 + MW2) );
    double jacobian = 1; //s * x; // from d2s/dxdQ2 --> d2s/dxdy    

    double F1_val = SF_PDF->xfxQ2(F1_code, x, Q2/SQ(pc->GeV));
    double F2_val = SF_PDF->xfxQ2(F2_code, x, Q2/SQ(pc->GeV));
    double F3_val = SF_PDF->xfxQ2(F3_code, x, Q2/SQ(pc->GeV));

    double F4_val = 0.0;
    double F5_val = F2_val / (2.0*x);
    
    double term1, term2, term3, term4, term5;
    if (config.XS.enable_mass_terms) {
        term1 = y * ( x*y + SQ(M_l)/(2*ENU*M_target) ) * F1_val;
        term2 = ( 1 - y - M_target*x*y/(2*ENU) - SQ(M_l)/(4*SQ(ENU)) ) * F2_val;
        term3 = config.cp_factor * (x*y*(1-y/2) - y*SQ(M_l)/(4*M_target*ENU)) * F3_val;
        term4 = (x*y*SQ(M_l) / (2 * M_target * ENU) + SQ(M_l)*SQ(M_l) / (4*SQ(M_target*ENU))) * F4_val;
        term5 = -1.0 * SQ(M_l) / (M_target * ENU) * F5_val;
    } else {
        term1 = y*y*x * F1_val;
        term2 = ( 1 - y ) * F2_val;
        term3 = config.cp_factor * ( x*y*(1-y/2) ) * F3_val;
        term4 = 0.0;
        term5 = 0.0;
    }
    double xs = fmax(prefactor * jacobian * propagator * (term1 + term2 + term3 + term4 + term5), 0);
    return xs / SQ(pc->cm); // TODO: Unit conversion outside of this function?
}

double CrossSection::ds_dxdy(double x, double y) {
    double xy = x*y;
    double MW2 = config.constants.Mboson2 * SQ(pc->GeV); // TODO: This should happen where M_boson2 is?
    double M_target = config.target_mass;
    double ME = M_target * ENU;
    double M_l2 = SQ(config.lepton_mass);

    double s = 2.0 * M_target * ENU + SQ(M_target);// - SQ(top_mass*pc->GeV);
    double Q2 = (s - SQ(M_target)) * xy;

    double prefactor = SQ(pc->GF) / (2 * M_PI * x); 
    double propagator = SQ( MW2 / (Q2 + MW2) );
    double jacobian = (s - SQ(M_target)) * x; // from d2s/dxdQ2 --> d2s/dxdy    

    double F1_val = SF_PDF->xfxQ2(F1_code, x, Q2/SQ(pc->GeV));
    double F2_val = SF_PDF->xfxQ2(F2_code, x, Q2/SQ(pc->GeV));
    double F3_val = SF_PDF->xfxQ2(F3_code, x, Q2/SQ(pc->GeV));

    // double F4_val = 0.0;
    double F5_val = F2_val / (2.0*x);
    
    // double term1, term2, term3, term4, term5;
    double term1, term2, term3, term5;
    if (config.XS.enable_mass_terms) {
        term1 = y * ( xy + M_l2/(2*ME) ) * F1_val;
        term2 = ( 1 - y - M_target*xy/(2*ENU) - M_l2/(4*SQ(ENU)) ) * F2_val;
        term3 = config.cp_factor * (xy*(1-y/2) - y*M_l2/(4*ME)) * F3_val;
        // term4 = (xy*M_l2 / (2 * ME) + M_l2*M_l2 / (4*SQ(ME))) * F4_val;
        term5 = -1.0 * M_l2 / (ME) * F5_val;
    } else {
        term1 = y*y*x * F1_val;
        term2 = ( 1 - y ) * F2_val;
        term3 = config.cp_factor * y*(1-y/2) * x*F3_val;
        // term4 = 0.0;
        term5 = 0.0;
    }
    // double xs = fmax(prefactor * jacobian * propagator * (term1 + term2 + term3 + term4 + term5), 0);
    double xs = fmax(prefactor * jacobian * propagator * (term1 + term2 + term3 + term5), 0);
    return xs / SQ(pc->cm); // TODO: Unit conversion outside of this function?
}

double CrossSection::ds_dy_kernel(double k) {
    double x = std::exp(k);

    bool ps_valid = ps.Validate(ENU, x, kernel_y);
    if (!ps_valid) {
        return 0;
    }

    double ddxs = ds_dxdy(x, kernel_y);
    if (config.XS.enable_radiative_corrections) {
        double radiative_correction = rc_dsdxdy(ENU, x, kernel_y, ddxs);
        ddxs += radiative_correction;
        // double radiative_correction = rc_bardin(ENU, x, kernel_y);
        // ddxs = ddxs * (1.0 + radiative_correction);
    }

    double result = x * ddxs;
    return result;
}


double CrossSection::ds_dy(double E, double y) {
    Set_Neutrino_Energy(E);
    double M_target = config.target_mass;
    double M_l = config.lepton_mass;

    kernel_y = y;
    // double s = 2.0 * M_target * E + SQ(M_target);
    double xmin = max(ps.x_min, SQ(M_l) / (2.0 * M_target * (E - M_l)));
    double xmax = ps.x_max;

    if (!ps.Validate(E)) {
        return 0;
    }

    gsl_integration_cquad_workspace * w = gsl_integration_cquad_workspace_alloc(2500);
    double result, error;
    size_t neval;

    gsl_function F;
    F.function = &KernelHelper<CrossSection, &CrossSection::ds_dy_kernel>;
    F.params = this;
    
    int status = gsl_integration_cquad(&F, log(xmin), log(xmax), 0, 1.e-4, w, &result, &error, &neval);
    if (status != 0) {
        std::cout << "ERR: " << status << std::endl;
    }
    gsl_integration_cquad_workspace_free(w);

    return result;
}

double CrossSection::TotalXS_xQ2(double E) {
    // integrate over x-Q2
    Set_Neutrino_Energy(E);
    double M_target = config.target_mass;
    double M_l = config.lepton_mass;

    double s = 2.0 * M_target * E + SQ(M_target);
    
    double xmin = max(ps.x_min, SQ(M_l) / (2.0 * M_target * (E - M_l)));
    double xmax = 1.0-1e-9;

    double Q2min = ps.Q2_min;
    double Q2max = s - SQ(M_target);

    double res,err;
    const unsigned long dim = 2; int calls = 250000; // bump it

    // integrating on the log of x and Q2
    double xl[dim] = { log(xmin), log(Q2min) };
    double xu[dim] = { log(xmax), log(Q2max) };

    gsl_rng_env_setup ();
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *r = gsl_rng_alloc (T);

    gsl_monte_function F;
    F = { &KernelHelper<CrossSection, &CrossSection::ds_dxdQ2_kernel>, dim, this};

    gsl_monte_vegas_state *s_vegas = gsl_monte_vegas_alloc (dim);
    gsl_monte_vegas_integrate (&F, xl, xu, dim, 10000, r, s_vegas, &res, &err);

    do
    {
        gsl_monte_vegas_integrate (&F, xl, xu, dim, calls/5, r, s_vegas, &res, &err);
        // printf ("result = % .6e sigma = % .6e "
        //         "chisq/dof = %.2f\n", res, err, gsl_monte_vegas_chisq (s_vegas));
    }
    while (fabs (gsl_monte_vegas_chisq (s_vegas) - 1.0) > 0.5 );

    gsl_monte_vegas_free (s_vegas);
    gsl_rng_free (r);

    return res;
}

double CrossSection::TotalXS(double E){
    Set_Neutrino_Energy(E);
    double M_target = config.target_mass;
    double M_l = config.lepton_mass;

    // double s = 2.0 * M_target * E + SQ(M_target);
    
    double xmin = max(ps.x_min, SQ(M_l) / (2.0 * M_target * (E - M_l)));
    double xmax = ps.x_max;

    if (!ps.Validate(E)) {
        return 0;
    }

    double res,err;
    const unsigned long dim = 2; int calls = 25000; // bump it

    // integrating on the log of x and y
    double xl[dim] = { log(xmin), log(1.e-12) };
    double xu[dim] = { log(xmax), log(1.)    };

    gsl_rng_env_setup ();
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *r = gsl_rng_alloc (T);

    gsl_monte_function F;
    F = { &KernelHelper<CrossSection, &CrossSection::ds_dxdy_kernel>, dim, this};
    gsl_monte_vegas_state *s_vegas = gsl_monte_vegas_alloc (dim);
    gsl_monte_vegas_integrate (&F, xl, xu, dim, 1000, r, s_vegas, &res, &err);

    int max_iterations = 10;
    int n_iter = 0;
    do
    {
        gsl_monte_vegas_integrate (&F, xl, xu, dim, calls/5, r, s_vegas, &res, &err);
        // printf ("result = % .6e sigma = % .6e "
        //         "chisq/dof = %.2f\n", res, err, gsl_monte_vegas_chisq (s_vegas));
        n_iter += 1;
    }
    while ( (fabs(gsl_monte_vegas_chisq (s_vegas) - 1.0) > 0.5) && (n_iter < max_iterations) );
    if (n_iter == max_iterations) {
        std::cout << "WARNING: Max iterations reached for E = " << E/pc->GeV << " GeV!" << std::endl;
    }

    gsl_monte_vegas_free (s_vegas);
    gsl_rng_free (r);

    return res;
}

/*
Radiative Corrections
*/
long double CrossSection::qed_splitting(long double z) {
    return (1.0L + z*z) / (1.0L - z);
    // return exp( log1p(z*z) - log1p(-z));
}

long double CrossSection::P11(long double z) {
    return qed_splitting(z);
}

long double CrossSection::P21(long double z) {
    return ((1.0L + z*z) / (1.0L - z)) * (2.0L * log(1-z) - log(z) + 1.5L) + 0.5L * (1.0L + z) * log(z) - (1.0 - z);
}

long double CrossSection::P22(long double z) {
    return (1.0L + z) * log(z) + 0.5L * (1.0L - z) + (2.0L / 3.0L) * (1.0L / z - z*z);
}

long double CrossSection::P23(long double z) {
    // leptons
}

long double CrossSection::Psoft(long double z) {
    long double eta = -3.0L * log(1.0L - (rc_prefactor / M_PI) * kernel_L );
    long double D = eta * pow(1.0L - z, eta - 1.0L) * exp(0.5L * eta * (1.5L - 2.0L * pc->gammaE)) / tgamma(1.0L+eta);
    return D - rc_prefactor * kernel_L * (2.0L / (1.0L - z)) * (1.0L + rc_prefactor * kernel_L * (11.0L/6.0L + 2.0L * log(1.0L-z)));
}

long double CrossSection::qed_splitting_integrated(long double a) {
    return -0.5L * a * (a + 2.0L) - 2.0L * std::log(1.0L - a);
}

long double CrossSection::rc_jacobian(long double x, long double y, long double z) {
     return y / (z * (z + y - 1.0L));
    // return exp( log(y) - log(z) - log(z + y - 1.0L));
}
double CrossSection::rc_kernel(double k) {
    return static_cast<double>(calculate_highorder_rc_dsdzdxdy(static_cast<long double>(k), 
                                                     static_cast<long double>(kernel_x), 
                                                     static_cast<long double>(kernel_y))
                              );
}

long double CrossSection::calculate_rc_dsdzdxdy(long double z, long double x, long double y) {
    long double xhat = x * y / (z + y - 1.0L);
    long double yhat = (z + y - 1.0L) / z;
    long double zmin = 1.0L - y * (1.0L - x);

    long double term1 = 0.0;
    long double born_dsdxdy_hat = 0.0;

    if (xhat > 1.0L) xhat = 1.0L;
    if (yhat > 1.0L) yhat = 1.0L;
    
    if (z > zmin) {
        if (ps.Validate(ENU, xhat, yhat)) { // only calc if in phase space, is this right? should entire dsdzdxdy be zero? -PW
            if (interp_grid_loaded) {
                double Ei[] = {ENU};
                double yi[] = {static_cast<double>(yhat)};
                double xi[] = {static_cast<double>(xhat)};
                double interp_output[1];
                mlinterp::interp(interp_nd, 1, interp_indata, interp_output, interp_E, Ei, interp_y, yi, interp_x, xi);
                born_dsdxdy_hat = static_cast<long double>(interp_output[0]);
            } else { // if we have no spline, do it manually
                born_dsdxdy_hat = static_cast<long double>(ds_dxdy(xhat, yhat));
            }
            
            if (born_dsdxdy_hat < 0) {
                born_dsdxdy_hat = 0.0L;
            }

            term1 = rc_jacobian(x, y, z) * born_dsdxdy_hat;
        }
    }

    double Ei[] = {ENU};
    double yi[] = {static_cast<double>(y)};
    double xi[] = {static_cast<double>(x)};
    double interp_output[1];
    mlinterp::interp(interp_nd, 1, interp_indata, interp_output, interp_E, Ei, interp_y, yi, interp_x, xi);
    long double born_dsdxdy = static_cast<long double>(interp_output[0]);

    // std::cout << "RC PS Failed: (E, x, y, z) = ( " << ENU/1e9 << ", " << xhat/x << ", " << yhat/y << ", " << zmin/z << " )" <<std::endl;
// std::cout << "            (xh, yh, zmin) = ( " << xhat << ", " << yhat << ", " << zmin << " )" <<std::endl;

    return qed_splitting(z) * (term1 - born_dsdxdy);
}

long double CrossSection::calculate_highorder_rc_dsdzdxdy(long double z, long double x, long double y) {
    long double xhat = x * y / (z + y - 1.0L);
    long double yhat = (z + y - 1.0L) / z;
    long double zmin = 1.0L - y * (1.0L - x);

    long double term1 = 0.0;
    long double born_dsdxdy_hat = 0.0;

    long double L = static_cast<long double>(kernel_L);
    long double a = static_cast<long double>(rc_prefactor);

    long double term_one = 0.0L;
    // long double term_two = 0.0L;
    // long double term_three = 0.0L;

    long double J = rc_jacobian(x, y, z);

    if (xhat > 1.0L) xhat = 1.0L;
    if (yhat > 1.0L) yhat = 1.0L;

    if (z > zmin) {
        
        if (ps.Validate(ENU/z, xhat, yhat)) { 
            // double Ei[] = {ENU/z};
            // double yi[] = {static_cast<double>(yhat)};
            // double xi[] = {static_cast<double>(xhat)};
            // double interp_output[1];
            // mlinterp::interp(interp_nd, 1, interp_indata, interp_output, interp_E, Ei, interp_y, yi, interp_x, xi);
            // born_dsdxdy_hat = static_cast<long double>(interp_output[0]);
            Set_Neutrino_Energy(ENU/z);
            born_dsdxdy_hat = static_cast<long double>(ds_dxdy(xhat, yhat));
            Set_Neutrino_Energy(ENU);
            term_one = J * born_dsdxdy_hat;
            // term_two = P21(z) * J * born_dsdxdy_hat;
            // term_three = P22(z) * J * born_dsdxdy_hat;
        }
    }

    // double Ei[] = {ENU};
    // double yi[] = {static_cast<double>(y)};
    // double xi[] = {static_cast<double>(x)};
    // double interp_output[1];
    // mlinterp::interp(interp_nd, 1, interp_indata, interp_output, interp_E, Ei, interp_y, yi, interp_x, xi);
    // long double born_dsdxdy = static_cast<long double>(interp_output[0]);
    long double born_dsdxdy = static_cast<long double>(kernel_born_dsdxdy);

    // return a*L*P11(z)*(term_one - born_dsdxdy) + a*a*L*L*(P21(z)*(term_one - born_dsdxdy) + P22(z)*term_one) + Psoft(z) * (term_one - born_dsdxdy);
    // return a*L*P11(z)*(term_one - born_dsdxdy) + a*a*L*L*(P21(z)*(term_one - born_dsdxdy) + P22(z)*term_one);
    return a*L*P11(z)*(term_one - born_dsdxdy);
}


double CrossSection::rc_dsdxdy(double E, double x, double y, double born_dsdxdy) {
    if (born_dsdxdy <= 0.0) return 0.0;

    Set_Neutrino_Energy(E);
    // TODO: think about if we want this to be called externally at all
    // In the kernel for dsdxdy, we call this function so ENU is already set and the phase space is checked
    double M_target = config.target_mass;
    double M_l = config.lepton_mass;
    kernel_y = y;
    kernel_x = x;
    kernel_born_dsdxdy = born_dsdxdy; // set this to speed up the calculation
    
    double s = 2.0 * M_target * E + SQ(M_target);
    double Q2 = (s - SQ(M_target)) * x * y;
    // kernel_L = log( (s * SQ(1.0 - y + x*y)) / SQ(M_l)); // large logarithm
    kernel_L = log( Q2 / SQ(M_l));
    // long double zmin = 1.0L - static_cast<long double>(y) * (1.0L - static_cast<long double>(x));

    double integrate_zmin = 1e-5;
    double integrate_zmax = 0.9999; // maybe this needs to be 0.9999 or something -PW
    if (!ps.Validate(E, x, y)) {
        return 0;
    }

    gsl_integration_cquad_workspace * w = gsl_integration_cquad_workspace_alloc(2500);
    double result, error;
    size_t neval;
    gsl_function F;
    F.function = &KernelHelper<CrossSection, &CrossSection::rc_kernel>;
    F.params = this;
    int status = gsl_integration_cquad(&F, integrate_zmin, integrate_zmax, 0, 1.e-3, w, &result, &error, &neval);
    if (status != 0) {
        std::cout << "RC ERR: (" << x << ", " << y << "): " << status << ", " << neval << " | " << born_dsdxdy << ", " << result << std::endl;
        result = 0.0;
    }
    gsl_integration_cquad_workspace_free(w);

    // if (abs(result) > 1e-60) {
    //     std::cout << result << std::endl;
    // }
    
    // return rc_prefactor * kernel_L * result;
    return result;
}

double CrossSection::rc_bardin(double E, double x, double y) {
    double M_target = config.target_mass;
    double M_l = config.lepton_mass;
    double MZ = pc->Zboson_mass;

    double s = 2.0 * M_target * E + SQ(M_target);
    double L = log(s / SQ(M_l));

    double zeta2 = 2.0*SQ(M_PI)/6.0;

    double delta_1LL = rc_prefactor * L * (2.0 * log(y) - log(1 - y) + 1.5 - y);
    double delta_2LL = 0.5 * SQ(rc_prefactor * L) * (4.0 * SQ(log(y)) - 4.0 * log(y) * log(1-y) + 0.5 * SQ(log(1-y) - 
                     4.0*zeta2 + (6.0 - 4.0*y)*log(y) + (3.0*y - 4.0)*log(1-y) + 9.0/4.0)) + (1.0/3.0)*rc_prefactor*L*delta_1LL;

    double prefactor = pc->alpha / M_PI;
    double line1 = -1.5 * log(s/SQ(MZ)) + (0.75 - 0.5*y - 0.5*log(1-y)+log(y)) * log(s/SQ(M_l));
    double line2 = 0.5 * log(1-y)*log(y) - 1.5*gsl_sf_dilog(y) - 0.5*SQ(log(y)) - 0.5*SQ(log(1-y));
    double line3 = (1-y)*log(1-y) - 7.0/4.0 * log(y) + 0.5*y*log(y) + 1.25*y + 0.5 + zeta2;

    double Q1 = -(2.0/3.0);
    double quark1 = Q1 * (3.0 - 1.5*log(s/SQ(MZ)) + zeta2 - log(1-y)*log(y) - 2*gsl_sf_dilog(y));
    double quark2 = SQ(Q1) * (23.0/72.0 -5.0/12.0 * y + SQ(y)/24.0 - zeta2 - 17.0/18.0);

    double delta_CC = prefactor * (line1 + line2 + line3 + quark1 + quark2);

    // if (delta_2LL != delta_2LL) {
    //     // std::cout << E << ", " << x << ", " << y << ": " << delta_1LL << ", " << delta_2LL << ", " << delta_CC << std::endl;
    //     std::cout << E << ", " << x << ", " << y << ": " << L << ", " << log(y) << ", " << log(1-y) << std::endl;
    // }


    return delta_CC + delta_2LL;
    // return delta_CC+delta_1LL;// + delta_2LL;
    // return delta_1LL + delta_2LL;
    // return delta_1LL;
}

void CrossSection::Load_RC_Spline(string spline_path) {
    if (rc_spline_loaded) {
        rc_spline = photospline::splinetable<>();
    }
    rc_spline.read_fits(spline_path);
    rc_spline_loaded = true;
}

void CrossSection::Load_InterpGrid(string grid_path) {

    std::cout << "Loading grid for linear interpolator: " << grid_path << std::endl;

    std::vector<double> E_values;
    std::vector<double> y_values;
    std::vector<double> x_values;

    std::ifstream gridfile(grid_path);
    std::string line;
    std::string discard;
    std::string token;

    // The first line contains the E values
    std::getline(gridfile, line);
    std::istringstream E_stream(line);
    while (std::getline(E_stream, token, ',')) {
        if (token == "E") continue;
        double E = std::stod(token);
        E_values.push_back(E);
    }

    // The second line contains the y values
    std::getline(gridfile, line);
    std::istringstream y_stream(line);
    while (std::getline(y_stream, token, ',')) {
        if (token == "y") continue;
        double y = std::stod(token);
        y_values.push_back(y);
    }

    // The second line contains the x values
    std::getline(gridfile, line);
    std::istringstream x_stream(line);
    while (std::getline(x_stream, token, ',')) {
        if (token == "x") continue;
        double x = std::stod(token);
        x_values.push_back(x);
    }

    int nE = E_values.size();
    int ny = y_values.size();
    int nx = x_values.size();

    int n = 0;
    for (int i = 0; i < nE; i++){
        for (int j = 0; j < ny; j++) {
            std::getline(gridfile, line);
            std::istringstream linestream(line);
            std::string val;
            int k = 0;
            while (std::getline(linestream, val, ',')) {
                n = k + j * nx + i * nx * ny;
                interp_indata[n] =  std::stod(val);
                interp_E[i] = E_values[i];
                interp_y[j] = y_values[j];
                interp_x[k] = x_values[k];
                k += 1;
            }
        }
    }
    gridfile.close();

    interp_grid_loaded = true;

    // int nE_interp = 1;
    // int nx_interp = 1;
    // int ny_interp = 1;

    // nd = {nE, ny, nx};
    // int ni = 1;

    // double E_interp[nE_interp] = {1e10};
    // double y_interp[ny_interp] = {0.1};
    // double x_interp[nx_interp] = {0.1};

    // double output[nE_interp*ny_interp*nx_interp]; // Result is stored in this buffer
    // mlinterp::interp(
    //     nd, ni,                   // Number of points
    //     data, output,                                          // Output axis
    //     data_E, E_interp, data_y, y_interp, data_x, x_interp);

    // for (auto& a : output) {
    //     std::cout << a << std::endl;
    // }
    // std::cout << nE << ", " << ny << ", " << nx << std::endl;

}

}