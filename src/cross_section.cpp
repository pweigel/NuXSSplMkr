#include "cross_section.h"

namespace nuxssplmkr {

CrossSection::CrossSection(Configuration& _config)  
    : config(_config)
{
    pc = new nuxssplmkr::PhysConst();
    M_iso = 0.5*(pc->proton_mass + pc->neutron_mass);

    // Set limits of integration
    // Note: these will change based on neutrino energy and certain features
    // integral_min_Q2 = config.Q2min;
    // integral_min_Q2 = 1.6;
    // integral_min_Q2 = 1.0;
    // if (config.enable_small_x) {
    //     integral_min_Q2 = 10.0;
    // }
    integral_min_x = config.XS.xmin;
    integral_max_x = config.XS.xmax;
    integral_min_Q2 = config.XS.Q2min;
    integral_max_Q2 = config.XS.Q2max;

}

void CrossSection::Load_Structure_Functions(string sf1_path, string sf2_path, string sf3_path) {
    Load_F1(sf1_path);
    Load_F2(sf2_path);
    Load_F3(sf3_path);
}

Grid CrossSection::Load_Grid(string path) {
    std::ifstream infile;
    infile.open(path);
    std::string line;

    Grid grid;
    std::vector<std::string> tokens;

    // Get first two lines
    std::getline(infile, line);
    boost::split(tokens, line, boost::is_any_of(" "));
    grid.NQ2 = stoi(tokens[0]);
    grid.Nx = stoi(tokens[1]);

    std::getline(infile, line);
    boost::split(tokens, line, boost::is_any_of(" "));
    grid.Q2min = stod(tokens[0]);
    grid.Q2max = stod(tokens[1]);
    grid.xmin = stod(tokens[2]);
    grid.xmax = stod(tokens[3]);

    int n = 0;
    double value;
    vector<double> data;
    while(std::getline(infile, line)) {
        std::stringstream linestream(line);
        std::string val;

        while(getline(linestream, val,',')) {
            data.push_back(stod(val));
        }
        n++;
    }
    grid.data = data;
    infile.close();

    return grid;
}

void CrossSection::Load_F1(string path) {
    bool is_grid = false;
    if (is_grid) {
        // Grid
        Grid grid = Load_Grid(path);
        double dlogQ2 = (grid.Q2max - grid.Q2min) / ((double)grid.NQ2 - 1);
        double dlogx = (grid.xmax - grid.xmin) / ((double)grid.Nx - 1);
        std::cout << "F1 Grid Limits: " << std::endl;
        std::cout << "    logQ2 = [" << grid.Q2min << ", " << grid.Q2max << "]" << std::endl;
        std::cout << "    logx = [" << grid.xmin << ", " << grid.xmax << "]" << std::endl;
        F1_interpolator =  BilinearInterpolator(std::move(grid.data), grid.NQ2, grid.Nx, grid.Q2min, grid.Q2max, grid.xmin, grid.xmax);

    } else {
        if (F1_loaded) {
            F1 = photospline::splinetable<>();
        }
        F1.read_fits(path);
    }
    std::cout << "Loaded F1!" << std::endl;
    F1_loaded = true;
}

void CrossSection::Load_F2(string path) {
    bool is_grid = false;
    if (is_grid) {
        // Grid
        Grid grid = Load_Grid(path);
        double dlogQ2 = (grid.Q2max - grid.Q2min) / ((double)grid.NQ2 - 1);
        double dlogx = (grid.xmax - grid.xmin) / ((double)grid.Nx - 1);
        F2_interpolator = BilinearInterpolator(std::move(grid.data), grid.NQ2, grid.Nx, grid.Q2min, grid.Q2max, grid.xmin, grid.xmax);

    } else {
        if (F2_loaded) {
            F2 = photospline::splinetable<>();
        }
        F2.read_fits(path);
    }
    std::cout << "Loaded F2!" << std::endl;
    F2_loaded = true;
}

void CrossSection::Load_F3(string path) {
    bool is_grid = false;
    if (is_grid) {
        // Grid
        Grid grid = Load_Grid(path);
        double dlogQ2 = (grid.Q2max - grid.Q2min) / ((double)grid.NQ2 - 1);
        double dlogx = (grid.xmax - grid.xmin) / ((double)grid.Nx - 1);
        F3_interpolator =  BilinearInterpolator(std::move(grid.data), grid.NQ2, grid.Nx, grid.Q2min, grid.Q2max, grid.xmin, grid.xmax);

    } else {
        if (F1_loaded) {
            F3 = photospline::splinetable<>();
        }
        F3.read_fits(path);
    }
    std::cout << "Loaded F3!" << std::endl;
    F3_loaded = true;
}

void CrossSection::Set_Neutrino_Energy(double E) {
    ENU = E;
    integral_max_Q2 = 2.0 * M_iso * ENU; // TODO: M_ISO --> target mass
}

void CrossSection::Set_Lepton_Mass(double m) {
    M_l = m;
}

double CrossSection::ds_dxdy_kernel(double* k) {
    double x = std::exp(k[0]);
    double y = std::exp(k[1]);

    // Integration limits
    double W2min = SQ(2.0 * pc->GeV); // TODO: This should be an input parameter
    double Q2 = 2.0 * M_iso * ENU * x * y;  // Q2 = s x y - m^2

    if ((x < 1e-9) || (x > 1)) {
        return 1e-99;
    }
    if (Q2 / SQ(pc->GeV) < integral_min_Q2) {
        return 1e-99;
    } 
    else if (Q2 / SQ(pc->GeV) > integral_max_Q2) {
        return 1e-99;
    }
    // double W2 = Q2 * (1.0 / x - 1.0) + SQ(M_iso); // TODO: target mass
    double W2 =  2.0 * M_iso * ENU * y * (1.0 - x) + SQ(M_iso); // Without the division // TODO: target mass

    // integration limits from Reno, Kretzer
    double _xmin = SQ(M_l) / (2 * M_iso * (ENU - M_l));

    double _a = (1.0 - SQ(M_l) * (1.0 / (2 * M_iso * ENU * x) + 1.0 / (2 * SQ(ENU))) ) / (2.0 * (1.0 + M_iso * x / (2.0 * ENU)));
    double _b = sqrt( SQ(1.0 - SQ(M_l)/(2.0 * M_iso * ENU * x)) ) / (2.0 * (1.0 + M_iso * x / (2.0 * ENU)));
    double _ymin = _a - _b;
    double _ymax = _a + _b;
    // std::cout << "ymin: " << _ymin << "," << "ymax: " << _ymax << std::endl;
    // if ((y < _ymin) || (y > _ymax)){
    //     return 1e-99;
    // }

    // TODO: better implementation
    if (W2 < W2min) {
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
    double MW2 = config.constants.Mboson2 * SQ(pc->GeV); // TODO: This should happen where M_boson2 is?

    double s = 2.0 * M_iso * ENU + SQ(M_iso);
    double Q2 = (s - SQ(M_iso)) * x * y;
    // std::cout << s << " " << Q2 << " " << x << " " << y;

    double prefactor = SQ(pc->GF) / (2 * M_PI * x); 
    double propagator = SQ( MW2 / (Q2 + MW2) );
    double jacobian = s * x; // from d2s/dxdQ2 --> d2s/dxdy    
    // std::cout  << "log10Q2, log10x = " << std::log10(Q2 / SQ(pc->GeV)) << ", " << std::log10(x) << std::endl;
    std::array<double, 2> pt{{std::log10(Q2 / SQ(pc->GeV)), std::log10(x)}};

    std::array<int, 2> F1_splc;
    std::array<int, 2> F2_splc;
    std::array<int, 2> F3_splc;

    F1.searchcenters(pt.data(), F1_splc.data());
    F2.searchcenters(pt.data(), F2_splc.data());
    F3.searchcenters(pt.data(), F3_splc.data());

    double F1_val = F1.ndsplineeval(pt.data(), F1_splc.data(), 0);
    // std::cout << F1_val << ", ";
    double F2_val = F2.ndsplineeval(pt.data(), F2_splc.data(), 0);
    // std::cout << F2_val << ", ";
    double F3_val = F3.ndsplineeval(pt.data(), F3_splc.data(), 0);
    // std::cout << F3_val << std::endl;
    
    double term1, term2, term3;
    if (config.XS.enable_mass_terms) {
        term1 = y*y*x + SQ(M_l) * y / (2 * ENU * M_iso) * F1_val;
        term2 = ((1.0 - SQ(M_l/ENU)/4.0) - (1.0 + M_iso * x / (2.0 * ENU)*y)) * F2_val;
        term3 = config.cp_factor * ( x*y*(1-y/2) - SQ(M_l) * y / (4.0 * ENU * M_iso) ) * F3_val;
    } else {
        term1 = y*y*x * F1_val;
        term2 = ( 1 - y ) * F2_val;
        term3 = config.cp_factor * ( x*y*(1-y/2) ) * F3_val;
    }
    double xs = prefactor * jacobian * propagator * (term1 + term2 + term3);
    return xs / SQ(pc->cm); // TODO: Unit conversion outside of this function?
}

// double CrossSection::ds_dxdy(double x, double y) {
//     double MW2 = config.constants.Mboson2 * SQ(pc->GeV); // TODO: This should happen where M_boson2 is?

//     double s = 2.0 * M_iso * ENU + SQ(M_iso);
//     double Q2 = (s - SQ(M_iso)) * x * y;

//     double prefactor = SQ(pc->GF) / (2 * M_PI * x); 
//     double propagator = SQ( MW2 / (Q2 + MW2) );
//     double jacobian = s * x; // from d2s/dxdQ2 --> d2s/dxdy    

//     double logQ2 = std::log10(Q2 / SQ(pc->GeV));
//     double logx = std::log10(x);
    
//     double F1_val = F1_interpolator.interpolate(logQ2, logx);
//     double F2_val = F2_interpolator.interpolate(logQ2, logx);
//     double F3_val = F3_interpolator.interpolate(logQ2, logx);

//     double term1, term2, term3;
//     if (config.XS.enable_mass_terms) {
//         term1 = y*y*x + SQ(M_l) * y / (2 * ENU * M_iso) * F1_val;
//         term2 = ((1.0 - SQ(M_l/ENU)/4.0) - (1.0 + M_iso * x / (2.0 * ENU)*y)) * F2_val;
//         term3 = config.cp_factor * ( x*y*(1-y/2) - SQ(M_l) * y / (4.0 * ENU * M_iso) ) * F3_val;
//     } else {
//         term1 = y*y*x * F1_val;
//         term2 = ( 1 - y ) * F2_val;
//         term3 = config.cp_factor * ( x*y*(1-y/2) ) * F3_val;
//     }
//     double xs = prefactor * jacobian * propagator * (term1 + term2 + term3);
//     return xs / SQ(pc->cm); // TODO: Unit conversion outside of this function?
// }

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
    
    double MW2 = config.constants.Mboson2 * SQ(pc->GeV); // TODO: This should happen where M_boson2 is?

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
    if (config.target == "proton") {
        d    = config.pdf.pdf -> xfxQ2(1, x, mQ2);
        dbar = config.pdf.pdf -> xfxQ2(-1, x, mQ2);
        u    = config.pdf.pdf -> xfxQ2(2, x, mQ2);
        ubar = config.pdf.pdf -> xfxQ2(-2, x, mQ2);
    } else if (config.target == "neutron") {
        d    = config.pdf.pdf -> xfxQ2(2, x, mQ2);
        dbar = config.pdf.pdf -> xfxQ2(-2, x, mQ2);
        u    = config.pdf.pdf -> xfxQ2(1, x, mQ2);
        ubar = config.pdf.pdf -> xfxQ2(-1, x, mQ2);
    } else {
        throw std::runtime_error("Unrecognized target type!");
    }

    double s    = config.pdf.pdf -> xfxQ2(3, x, mQ2);
    double sbar = config.pdf.pdf -> xfxQ2(-3, x, mQ2);
    double b    = config.pdf.pdf -> xfxQ2(5, x, mQ2);
    double bbar = config.pdf.pdf -> xfxQ2(-5, x, mQ2);

    double c    = config.pdf.pdf -> xfxQ2(4, x, mQ2);
    double cbar = config.pdf.pdf -> xfxQ2(-4, x, mQ2);
    double t    = config.pdf.pdf -> xfxQ2(6, x, mQ2);
    double tbar = config.pdf.pdf -> xfxQ2(-6, x, mQ2);

    double F2_val;
    double xF3_val; // TODO:
    if (config.sf_type == SFType::charm) {
        if (config.neutrino_type == NeutrinoType::neutrino) {
            F2_val =  2 * (SQ(config.constants.Vcd)*d + SQ(config.constants.Vcs)*s + SQ(config.constants.Vcb)*b);
            xF3_val = 2 * (SQ(config.constants.Vcd)*d + SQ(config.constants.Vcs)*s + SQ(config.constants.Vcb)*b);
        } else if (config.neutrino_type == NeutrinoType::antineutrino) {
            F2_val =   2 * (SQ(config.constants.Vcd)*dbar + SQ(config.constants.Vcs)*sbar + SQ(config.constants.Vcb)*bbar);
            xF3_val = -2 * (SQ(config.constants.Vcd)*dbar + SQ(config.constants.Vcs)*sbar + SQ(config.constants.Vcb)*bbar);
        } else {
            throw std::runtime_error("Unrecognized neutrino type!");
        }
    } else {
        if (config.neutrino_type == NeutrinoType::neutrino) {
            F2_val =  2 * (d + s + b + ubar + cbar + tbar);
            xF3_val = 2 * (d + s + b - ubar - cbar - tbar);
        } else if (config.neutrino_type == NeutrinoType::antineutrino) {
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
    double term3 = config.cp_factor * ( y*(1-y/2) ) * xF3_val; // note: removed x here because we are computing xF3
    double xs = prefactor * jacobian * propagator * (term1 + term2 + term3);
    return xs / SQ(pc->cm); // TODO: Unit conversion outside of this function?
}

double CrossSection::_ds_dy(double k) {
    double x = std::exp(k);
    double Q2 = 2.0 * M_iso * ENU * x * _kernel_y;  // Q2 = s x y - m^2
    if ((x < config.SF.xmin) || (x > config.SF.xmax)) {
        return 1e-99;
    }
    if (Q2 / SQ(pc->GeV) < 1.5) {
        return 1e-99;
    }
    // std::cout << x << ", " << Q2 / SQ(pc->GeV) << std::endl;
    // if ((Q2 / SQ(pc->GeV) < config.SF.Q2min) || (Q2 / SQ(pc->GeV) > config.SF.Q2max)) {
    //     return 1e-99;
    // }

    // double _xmin = SQ(M_l) / (2 * M_iso * (ENU - M_l));
    // std::cout << "XMIN = " << _xmin << std::endl;

    // Integration limits
    // double _a = (1.0 - SQ(M_l) * (1.0 / (2 * M_iso * ENU * x) + 1.0 / (2 * SQ(ENU))) ) / (2.0 * (1.0 + M_iso * x / (2.0 * ENU)));
    // double _b = sqrt( SQ(1.0 - SQ(M_l)/(2 * M_iso * ENU * x)) ) / (2.0 * (1.0 + M_iso * x / (2.0 * ENU)));
    
    // double _ymin = _a - _b;
    // double _ymax = _a + _b;
    // std::cout << "YMIN = " << _ymin << ", " << "YMAX = " << _ymax << std::endl;

    double W2min;
    if (config.sf_type == charm) {
        W2min = SQ( (0.938 + 1.869) * pc->GeV); // (m_N + m_D)^2
    } else {
        W2min = SQ(2.0 * pc->GeV); // TODO: This should be an input parameter
    }
    // double W2 = Q2 * (1.0 / x - 1.0) + SQ(M_iso); // TODO: target mass
    // double W2 =  2.0 * M_iso * ENU * _kernel_y * (1.0 - x) + SQ(M_iso); // Without the division // TODO: target mass
    // std::cout << x << " " << _kernel_y << std::endl;
    // std::cout << "Q2: " << Q2 / SQ(pc->GeV) << " MINIMUM: " << integral_min_Q2 << std::endl;
    // std::cout << "W2: " << W2 / SQ(pc->GeV) << " MINIMUM: " << W2min / SQ(pc->GeV) << std::endl;
    // TODO: better implementation
    // if ((_kernel_y < _ymin) || (_kernel_y > _ymax)){
    //     return 1e-99;
    // }
    // if (W2 < W2min) {
    //     return 1e-99;
    // }
    // if (Q2 / SQ(pc->GeV) < integral_min_Q2) {
    //     return 1e-99;
    // } 
    // else if (Q2 / SQ(pc->GeV) > integral_max_Q2) {
    //     return 1e-99;
    // }
    
    double result = x * ds_dxdy(x, _kernel_y);
    // std::cout << x << ", " << Q2 << " " << result << std::endl;
    return result;
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

    double W2min = SQ(2.0 * pc->GeV);
    // double Q2 = 2.0 * M_iso * E * x * _kernel_y;
    double xmax = SQ(pc->GeV) * integral_max_Q2 / (2.0 * M_iso * E * y);
    double xmin = SQ(pc->GeV) * integral_min_Q2 / (2.0 * M_iso * E * y);
    if (xmin > 1) {
        return 0.0;
    }
    if (xmax > 1) {
        xmax = 1.0;
    }

    // double Q2min = 2.0 * M_iso * E * xmin * _kernel_y;
    // q2min = 2 gev
    // xmin = 1/m
    // double Q2max = 2.0 * M_iso * E * xmax * _kernel_y;
    // if (Q2min  / SQ(pc->GeV) < 2.0) {
    //     return 0.0;
    // }
    // std::cout << xmin << " " << xmax << std::endl;
    // std::cout << Q2min / SQ(pc->GeV) << " " << Q2max / SQ(pc->GeV) << std::endl;

    // double W2 =  2.0 * M_iso * E * _kernel_y * (1.0 - x) + SQ(M_iso); // Without the division // TODO: target mass
    // if (W2 < W2min) {
    //     return 1e-99;
    // }
    // if (Q2 / SQ(pc->GeV) < integral_min_Q2) {
    //     return 1e-99;
    // } 
    // else if (Q2 / SQ(pc->GeV) > integral_max_Q2) {
    //     return 1e-99;
    // }
    

    gsl_integration_workspace * w = gsl_integration_workspace_alloc(10000);
    double result, error;

    gsl_function F;
    if (config.SF.mass_scheme == "parton") { 
        F.function = &KernelHelper<CrossSection, &CrossSection::_ds_dy_partonic>;
    } else {
        F.function = &KernelHelper<CrossSection, &CrossSection::_ds_dy>;
    }
    F.params = this;
    
    gsl_integration_qag ( &F, log(xmin), log(xmax), 0, 1.e-5, 10000, 6, w, &result, &error);
    gsl_integration_workspace_free(w);

    return result;
}

// double CrossSection::ds_dy(double E, double y) {
//     Set_Neutrino_Energy(E);
//     _kernel_y = y

//     double res,err;
//     const unsigned long dim = 1; int calls = 50000;

//     // integrating on the log of x and y
//     double xl[dim] = { log(integral_min_x) };
//     double xu[dim] = { log(integral_max_x) };

//     gsl_rng_env_setup ();
//     const gsl_rng_type *T = gsl_rng_default;
//     gsl_rng *r = gsl_rng_alloc (T);

//     gsl_monte_function F;
//     F.f = &KernelHelper<CrossSection, &CrossSection::_ds_dy>;
//     F.dim = 1;
//     F.params = this;
//     // F = { &KernelHelper<CrossSection, &CrossSection::_ds_dy>, dim, this};

//     gsl_monte_vegas_state *s_vegas = gsl_monte_vegas_alloc (dim);

//     gsl_monte_vegas_integrate (&F, xl, xu, dim, 10000, r, s_vegas, &res, &err);

//     do
//     {
//         gsl_monte_vegas_integrate (&F, xl, xu, dim, calls, r, s_vegas, &res, &err);
//     }
//     while (fabs (gsl_monte_vegas_chisq (s_vegas) - 1.0) > 0.5 );

//     gsl_monte_vegas_free (s_vegas);
//     gsl_rng_free (r);

//     return res;
// }

double CrossSection::TotalXS(double E){
    Set_Neutrino_Energy(E);

    double res,err;
    const unsigned long dim = 2; int calls = 50000;

    // integrating on the log of x and y
    double xl[dim] = { log(integral_min_x), log(1.e-9) };
    double xu[dim] = { log(integral_max_x)   , log(1.)    };

    gsl_rng_env_setup ();
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *r = gsl_rng_alloc (T);

    gsl_monte_function F;
    if (config.SF.mass_scheme == "parton") {
        F = { &KernelHelper<CrossSection, &CrossSection::_ds_dxdy_partonic>, dim, this};
    } else {
        F = { &KernelHelper<CrossSection, &CrossSection::ds_dxdy_kernel>, dim, this};
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