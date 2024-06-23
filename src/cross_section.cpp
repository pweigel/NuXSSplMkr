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
    double F5_val = F2_val / x;
    
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
    double jacobian = s * x; // from d2s/dxdQ2 --> d2s/dxdy    

    double F1_val = SF_PDF->xfxQ2(F1_code, x, Q2/SQ(pc->GeV));
    double F2_val = SF_PDF->xfxQ2(F2_code, x, Q2/SQ(pc->GeV));
    double F3_val = SF_PDF->xfxQ2(F3_code, x, Q2/SQ(pc->GeV));

    // double F4_val = 0.0;
    double F5_val = F2_val / x;
    
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
    const unsigned long dim = 2; int calls = 250000; // bump it

    // integrating on the log of x and y
    double xl[dim] = { log(xmin), log(1.e-12) };
    double xu[dim] = { log(xmax), log(1.)    };

    gsl_rng_env_setup ();
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *r = gsl_rng_alloc (T);

    gsl_monte_function F;
    F = { &KernelHelper<CrossSection, &CrossSection::ds_dxdy_kernel>, dim, this};
    gsl_monte_vegas_state *s_vegas = gsl_monte_vegas_alloc (dim);
    gsl_monte_vegas_integrate (&F, xl, xu, dim, 10000, r, s_vegas, &res, &err);

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
double CrossSection::qed_splitting(double z) {
    return (1.0 + z*z) / (1.0 - z);
    // return exp( log1p(z*z) - log1p(-z) );
}

double CrossSection::rc_jacobian(double x, double y, double z) {
    return y / (z * (z + y - 1.0));
    // return exp( log(y) - log(z) - log(z + y - 1.0) );
}

double CrossSection::rc_kernel(double k) {
    return calculate_rc_dsdzdxdy(k, kernel_x, kernel_y);
}

// double CrossSection::calculate_rc_dsdxdy(double z, double E, double x, double y) {
//     double M_target = config.target_mass;
//     double M_l = config.lepton_mass;
//     double s = 2.0 * M_target * E + SQ(M_target);
//     double L = log( (s * SQ(1.0 - y + x*y)) / SQ(M_l));
//     return calculate_rc_dsdxdy(z, E, x, y, L);
// }

double CrossSection::calculate_rc_dsdzdxdy(double z, double x, double y) {
    double xhat = x * y / (z + y - 1.0);
    double yhat = (z + y - 1.0) / z;
    double zmin = 1.0 - y * (1.0 - x);

    // double M_target = config.target_mass;
    // double M_l = config.lepton_mass;
    // double s = 2.0 * M_target * E + SQ(M_target);
    // double L = log( (s * SQ(1.0 - y + x*y)) / SQ(M_l)); // large logarithm
    // double born_dsdxdy = ds_dxdy(x, y); // cross section evaluated at regular x, y
    double term1 = 0.0;
    if (z > zmin) {
        // if (xhat >= 1.0) xhat = 1.0;
        // if (yhat >= 1.0) yhat = 1.0;
        if (ps.Validate(ENU, xhat, yhat)) { // only calc if in phase space, is this right? should entire dsdzdxdy be zero? -PW
            double born_dsdxdy_hat = ds_dxdy(xhat, yhat);
            term1 = rc_jacobian(x, y, z) * born_dsdxdy_hat;
        }
    }

    return qed_splitting(z) * (term1 - kernel_born_dsdxdy);
}

double CrossSection::rc_dsdxdy(double E, double x, double y, double born_dsdxdy) {
    Set_Neutrino_Energy(E);
    // TODO: think about if we want this to be called externally at all
    // In the kernel for dsdxdy, we call this function so ENU is already set and the phase space is checked
    double M_target = config.target_mass;
    double M_l = config.lepton_mass;
    kernel_y = y;
    kernel_x = x;
    kernel_born_dsdxdy = born_dsdxdy; // set this to speed up the calculation
    
    double s = 2.0 * M_target * E + SQ(M_target);
    kernel_L = log( (s * SQ(1.0 - y + x*y)) / SQ(M_l)); // large logarithm

    double integrate_zmin = 0;
    double integrate_zmax = 0.99999; // maybe this needs to be 0.9999 or something -PW

    // E -> E/z > E

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
        std::cout << "RC ERR: " << status << std::endl;
    }
    gsl_integration_cquad_workspace_free(w);

    return rc_prefactor * kernel_L * result;
}

}