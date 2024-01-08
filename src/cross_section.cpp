#include "cross_section.h"

namespace nuxssplmkr {

CrossSection::CrossSection(Configuration& config) {
    sf_info = config.sf_info;
    pc = new nuxssplmkr::PhysConst();
    M_iso = 0.5*(pc->proton_mass + pc->neutron_mass);

    // Set limits of integration
    // Note: these will change based on neutrino energy and certain features
    // integral_min_Q2 = sf_info.Q2min;
    integral_min_Q2 = 1.6;
    if (sf_info.enable_small_x) {
        integral_min_Q2 = 10.0;
    }
    integral_max_Q2 = sf_info.Q2max;
    integral_min_x = sf_info.xmin;
    integral_max_x = sf_info.xmax;
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
    integral_max_Q2 = 2.0 * M_iso * ENU; // TODO: M_ISO --> target mass
}

double CrossSection::_ds_dxdy(double* k) {
    double x = std::exp(k[0]);
    double y = std::exp(k[1]);

    // Integration limits
    double W2min = SQ(2.0 * pc->GeV); // TODO: This should be an input parameter
    double Q2 = 2.0 * M_iso * ENU * x * y;  // Q2 = s x y - m^2
    // double W2 = Q2 * (1.0 / x - 1.0) + SQ(M_iso); // TODO: target mass
    double W2 =  2.0 * M_iso * ENU * y * (1.0 - x) + SQ(M_iso); // Without the division // TODO: target mass

    // TODO: better implementation
    if (W2 < W2min) {
        return 1e-99;
    }

    if (Q2 / SQ(pc->GeV) < integral_min_Q2) {
        return 1e-99;
    } 
    else if (Q2 / SQ(pc->GeV) > integral_max_Q2) {
        return 1e-99;
    }

    // Jacobian = x * y, needed because we're integrating over log space
    return x * y * ds_dxdy(x, y);
}

double CrossSection::ds_dxdy(double E, double x, double y) {
    Set_Neutrino_Energy(E);
    return ds_dxdy(x, y);
}

double CrossSection::ds_dxdy(double x, double y) {
    double MW2 = sf_info.M_boson2 * SQ(pc->GeV); // TODO: This should happen where M_boson2 is?

    double s = 2.0 * M_iso * ENU + SQ(M_iso); // TODO: target mass
    double Q2 = (s - SQ(M_iso)) * x * y; // TODO: target mass

    double prefactor = SQ(pc->GF) / (2 * M_PI * x); 
    double propagator = SQ( MW2 / (Q2 + MW2) );
    double jacobian = s * x; // from d2s/dxdQ2 --> d2s/dxdy    
    
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
    double term3 = sf_info.cp_factor * ( x*y*(1-y/2) ) * F3_val;

    double xs = prefactor * jacobian * propagator * (term1 + term2 + term3);
    return xs / SQ(pc->cm); // TODO: Unit conversion outside of this function?
}

double CrossSection::_ds_dxdy_partonic(double* k) {
    double x = std::exp(k[0]);
    double y = std::exp(k[1]);
    return x * y * ds_dxdy_partonic(x, y);
}

double CrossSection::ds_dxdy_partonic(double E, double x, double y) {
    Set_Neutrino_Energy(E);
    return ds_dxdy_partonic(x, y);
}

double CrossSection::ds_dxdy_partonic(double x, double y) {
    // TODO: W threshold and slow rescaling, CP factor fix
    
    double MW2 = sf_info.M_boson2 * SQ(pc->GeV); // TODO: This should happen where M_boson2 is?

    double s_energy = 2.*M_iso*ENU + SQ(M_iso); // using s for strange parton later
    double Q2 = (s_energy - SQ(M_iso)) * x * y;
    double mQ2 = Q2/(pc->GeV2); // This is used a lot, so precompute it

    double prefactor = SQ(pc->GF) / (2 * M_PI * x); 
    double propagator = SQ( MW2 / (Q2 + MW2) );
    double jacobian = s_energy * x; // from d2s/dxdQ2 --> d2s/dxdy

    // TODO: Generalize p/n and nu/nubar
    double d;
    double dbar;
    double u;
    double ubar;
    if (sf_info.target == "proton") {
        d    = sf_info.pdf -> xfxQ2(1, x, mQ2);
        dbar = sf_info.pdf -> xfxQ2(-1, x, mQ2);
        u    = sf_info.pdf -> xfxQ2(2, x, mQ2);
        ubar = sf_info.pdf -> xfxQ2(-2, x, mQ2);
    } else if (sf_info.target == "neutron") {
        d    = sf_info.pdf -> xfxQ2(2, x, mQ2);
        dbar = sf_info.pdf -> xfxQ2(-2, x, mQ2);
        u    = sf_info.pdf -> xfxQ2(1, x, mQ2);
        ubar = sf_info.pdf -> xfxQ2(-1, x, mQ2);
    } else {
        throw std::runtime_error("Unrecognized target type!");
    }

    double s    = sf_info.pdf -> xfxQ2(3, x, mQ2);
    double sbar = sf_info.pdf -> xfxQ2(-3, x, mQ2);
    double b    = sf_info.pdf -> xfxQ2(5, x, mQ2);
    double bbar = sf_info.pdf -> xfxQ2(-5, x, mQ2);

    double c    = sf_info.pdf -> xfxQ2(4, x, mQ2);
    double cbar = sf_info.pdf -> xfxQ2(-4, x, mQ2);
    double t    = sf_info.pdf -> xfxQ2(6, x, mQ2);
    double tbar = sf_info.pdf -> xfxQ2(-6, x, mQ2);

    double F2_val;
    double xF3_val; // TODO:
    if (sf_info.sf_type == SFType::charm) {
        if (sf_info.neutrino_type == NeutrinoType::neutrino) {
            F2_val =  2 * (SQ(sf_info.Vcd)*d + SQ(sf_info.Vcs)*s + SQ(sf_info.Vcb)*b);
            xF3_val = 2 * (SQ(sf_info.Vcd)*d + SQ(sf_info.Vcs)*s + SQ(sf_info.Vcb)*b);
        } else if (sf_info.neutrino_type == NeutrinoType::antineutrino) {
            F2_val =   2 * (SQ(sf_info.Vcd)*dbar + SQ(sf_info.Vcs)*sbar + SQ(sf_info.Vcb)*bbar);
            xF3_val = -2 * (SQ(sf_info.Vcd)*dbar + SQ(sf_info.Vcs)*sbar + SQ(sf_info.Vcb)*bbar);
        } else {
            throw std::runtime_error("Unrecognized neutrino type!");
        }
    } else {
        if (sf_info.neutrino_type == NeutrinoType::neutrino) {
            F2_val =  2 * (d + s + b + ubar + cbar + tbar);
            xF3_val = 2 * (d + s + b - ubar - cbar - tbar);
        } else if (sf_info.neutrino_type == NeutrinoType::antineutrino) {
            F2_val =  2 * (u + c + t + dbar + sbar + bbar);
            xF3_val = 2 * (u + c + t - dbar - sbar - bbar);
        } else {
            throw std::runtime_error("Unrecognized neutrino type!");
        }
    }

    double F1_val = F2_val / (2.0 * x);

    double term1 = y * ( y * x ) * F1_val;
    double term2 = ( 1 - y ) * F2_val;
    // double term3 = ( x*y*(1-y/2) ) * F3_val;
    double term3 = sf_info.cp_factor * ( y*(1-y/2) ) * xF3_val; // note: removed x here because we are computing xF3
    double xs = prefactor * jacobian * propagator * (term1 + term2 + term3);
    return xs / SQ(pc->cm); // TODO: Unit conversion outside of this function?
}

double CrossSection::_ds_dy(double k) {
    double x = std::exp(k);

    // Integration limits
    double W2min;
    if (sf_info.sf_type == charm) {
        W2min = SQ( (0.938 + 1.869) * pc->GeV); // (m_N + m_D)^2
    } else {
        W2min = SQ(2.0 * pc->GeV); // TODO: This should be an input parameter
    }
    double Q2 = 2.0 * M_iso * ENU * x * _kernel_y;  // Q2 = s x y - m^2
    // double W2 = Q2 * (1.0 / x - 1.0) + SQ(M_iso); // TODO: target mass
    double W2 =  2.0 * M_iso * ENU * _kernel_y * (1.0 - x) + SQ(M_iso); // Without the division // TODO: target mass

    // TODO: better implementation
    if (W2 < W2min) {
        return 1e-99;
    }

    if (Q2 / SQ(pc->GeV) < integral_min_Q2) {
        return 1e-99;
    } 
    else if (Q2 / SQ(pc->GeV) > integral_max_Q2) {
        return 1e-99;
    }
    return x * ds_dxdy(x, _kernel_y);
}

double CrossSection::_ds_dy_partonic(double k) {
    double x = std::exp(k);
    return x * ds_dxdy_partonic(x, _kernel_y);
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

    gsl_integration_workspace * w = gsl_integration_workspace_alloc(10000);
    double result, error;

    gsl_function F;
    if (sf_info.mass_scheme == "parton") { 
        F.function = &KernelHelper<CrossSection, &CrossSection::_ds_dy_partonic>;
    } else {
        F.function = &KernelHelper<CrossSection, &CrossSection::_ds_dy>;
    }
    F.params = this;
    
    gsl_integration_qag ( &F, log(1.e-9), log(1.), 0, 1.e-5, 10000, 6, w, &result, &error);
    gsl_integration_workspace_free(w);

    return result;
}

double CrossSection::TotalXS(double E){
    Set_Neutrino_Energy(E);

    double res,err;
    const unsigned long dim = 2; int calls = 50000;

    // integrating on the log of x and y
    double xl[dim] = { log(1.e-9), log(1.e-9) };
    double xu[dim] = { log(1.)   , log(1.)    };

    gsl_rng_env_setup ();
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *r = gsl_rng_alloc (T);

    gsl_monte_function F;
    if (sf_info.mass_scheme == "parton") {
        F = { &KernelHelper<CrossSection, &CrossSection::_ds_dxdy_partonic>, dim, this};
    } else {
        F = { &KernelHelper<CrossSection, &CrossSection::_ds_dxdy>, dim, this};
    }
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