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

// void CrossSection::Load_Structure_Functions(string sf1_path, string sf2_path, string sf3_path) {
//     Load_F1(sf1_path);
//     Load_F2(sf2_path);
//     Load_F3(sf3_path);
// }

// void CrossSection::Load_Structure_Functions(string inpath) {
//     string sf1_path = inpath + "/F1_"+config.projectile+"_"+config.target+"_"+config.sf_type_string+insuffix+".fits";
//     string sf2_path = inpath + "/F2_"+config.projectile+"_"+config.target+"_"+config.sf_type_string+insuffix+".fits";
//     string sf3_path = inpath + "/F3_"+config.projectile+"_"+config.target+"_"+config.sf_type_string+insuffix+".fits";
//     std::cout << sf1_path << std::endl;
//     std::cout << sf2_path << std::endl;
//     std::cout << sf3_path << std::endl;
//     Load_F1(sf1_path);
//     Load_F2(sf2_path);
//     Load_F3(sf3_path);
// }

// Grid CrossSection::Load_Grid(string path) {
//     std::ifstream infile;
//     infile.open(path);
//     std::string line;

//     Grid grid;
//     std::vector<std::string> tokens;

//     // Get first two lines
//     std::getline(infile, line);
//     boost::split(tokens, line, boost::is_any_of(" "));
//     grid.NQ2 = stoi(tokens[0]);
//     grid.Nx = stoi(tokens[1]);

//     std::getline(infile, line);
//     boost::split(tokens, line, boost::is_any_of(" "));
//     grid.Q2min = stod(tokens[0]);
//     grid.Q2max = stod(tokens[1]);
//     grid.xmin = stod(tokens[2]);
//     grid.xmax = stod(tokens[3]);

//     int n = 0;
//     double value;
//     vector<double> data;
//     while(std::getline(infile, line)) {
//         std::stringstream linestream(line);
//         std::string val;

//         while(getline(linestream, val,',')) {
//             data.push_back(stod(val));
//         }
//         n++;
//     }
//     grid.data = data;
//     infile.close();

//     return grid;
// }

// void CrossSection::Load_F1(string path) {
//     bool is_grid = false;
//     if (is_grid) {
//         // Grid
//         Grid grid = Load_Grid(path);
//         double dlogQ2 = (grid.Q2max - grid.Q2min) / ((double)grid.NQ2 - 1);
//         double dlogx = (grid.xmax - grid.xmin) / ((double)grid.Nx - 1);
//         std::cout << "F1 Grid Limits: " << std::endl;
//         std::cout << "    logQ2 = [" << grid.Q2min << ", " << grid.Q2max << "]" << std::endl;
//         std::cout << "    logx = [" << grid.xmin << ", " << grid.xmax << "]" << std::endl;
//         F1_interpolator =  BilinearInterpolator(std::move(grid.data), grid.NQ2, grid.Nx, grid.Q2min, grid.Q2max, grid.xmin, grid.xmax);

//     } else {
//         if (F1_loaded) {
//             F1 = photospline::splinetable<>();
//         }
//         F1.read_fits(path);
//     }
//     // std::cout << "Loaded F1!" << std::endl;
//     F1_loaded = true;
// }

// void CrossSection::Load_F2(string path) {
//     bool is_grid = false;
//     if (is_grid) {
//         // Grid
//         Grid grid = Load_Grid(path);
//         double dlogQ2 = (grid.Q2max - grid.Q2min) / ((double)grid.NQ2 - 1);
//         double dlogx = (grid.xmax - grid.xmin) / ((double)grid.Nx - 1);
//         F2_interpolator = BilinearInterpolator(std::move(grid.data), grid.NQ2, grid.Nx, grid.Q2min, grid.Q2max, grid.xmin, grid.xmax);

//     } else {
//         if (F2_loaded) {
//             F2 = photospline::splinetable<>();
//         }
//         F2.read_fits(path);
//     }
//     // std::cout << "Loaded F2!" << std::endl;
//     F2_loaded = true;
// }

// void CrossSection::Load_F3(string path) {
//     bool is_grid = false;
//     if (is_grid) {
//         // Grid
//         Grid grid = Load_Grid(path);
//         double dlogQ2 = (grid.Q2max - grid.Q2min) / ((double)grid.NQ2 - 1);
//         double dlogx = (grid.xmax - grid.xmin) / ((double)grid.Nx - 1);
//         F3_interpolator =  BilinearInterpolator(std::move(grid.data), grid.NQ2, grid.Nx, grid.Q2min, grid.Q2max, grid.xmin, grid.xmax);

//     } else {
//         if (F1_loaded) {
//             F3 = photospline::splinetable<>();
//         }
//         F3.read_fits(path);
//     }
//     // std::cout << "Loaded F3!" << std::endl;
//     F3_loaded = true;
// }

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

    // Jacobian = x * y, needed because we're integrating over log space
    return x * y * ds_dxdy(x, y);
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
    // std::array<double, 2> pt{{std::log10(Q2 / SQ(pc->GeV)), std::log10(x)}};

    // std::array<int, 2> F1_splc;
    // std::array<int, 2> F2_splc;
    // std::array<int, 2> F3_splc;

    // F1.searchcenters(pt.data(), F1_splc.data());
    // F2.searchcenters(pt.data(), F2_splc.data());
    // F3.searchcenters(pt.data(), F3_splc.data());

    // double F1_val = F1.ndsplineeval(pt.data(), F1_splc.data(), 0);
    // double F2_val = F2.ndsplineeval(pt.data(), F2_splc.data(), 0);
    // double F3_val = F3.ndsplineeval(pt.data(), F3_splc.data(), 0);

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
    double MW2 = config.constants.Mboson2 * SQ(pc->GeV); // TODO: This should happen where M_boson2 is?
    double M_target = config.target_mass;
    double M_l = config.lepton_mass;

    double s = 2.0 * M_target * ENU + SQ(M_target);// - SQ(top_mass*pc->GeV);
    double Q2 = (s - SQ(M_target)) * x * y;

    double prefactor = SQ(pc->GF) / (2 * M_PI * x); 
    double propagator = SQ( MW2 / (Q2 + MW2) );
    double jacobian = s * x; // from d2s/dxdQ2 --> d2s/dxdy    
    // std::array<double, 2> pt{{std::log10(Q2 / SQ(pc->GeV)), std::log10(x)}};

    // std::array<int, 2> F1_splc;
    // std::array<int, 2> F2_splc;
    // std::array<int, 2> F3_splc;

    // F1.searchcenters(pt.data(), F1_splc.data());
    // F2.searchcenters(pt.data(), F2_splc.data());
    // F3.searchcenters(pt.data(), F3_splc.data());

    // double F1_val = F1.ndsplineeval(pt.data(), F1_splc.data(), 0);
    // double F2_val = F2.ndsplineeval(pt.data(), F2_splc.data(), 0);
    // double F3_val = F3.ndsplineeval(pt.data(), F3_splc.data(), 0);

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
        term3 = config.cp_factor * y*(1-y/2) * x*F3_val;
        term4 = 0.0;
        term5 = 0.0;
    }
    double xs = fmax(prefactor * jacobian * propagator * (term1 + term2 + term3 + term4 + term5), 0);
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
    double MW2 = config.constants.Mboson2 * SQ(pc->GeV); // TODO: This should happen where M_boson2 is?
    double M_target = config.target_mass;

    double s_energy = 2.*M_target*ENU + SQ(M_target); // using s for strange parton later
    double Q2 = (s_energy - SQ(M_target)) * x * y;
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
        d    = config.pdf -> xfxQ2(1, x, mQ2);
        dbar = config.pdf -> xfxQ2(-1, x, mQ2);
        u    = config.pdf -> xfxQ2(2, x, mQ2);
        ubar = config.pdf -> xfxQ2(-2, x, mQ2);
    } else if (config.target == "neutron") {
        d    = config.pdf -> xfxQ2(2, x, mQ2);
        dbar = config.pdf -> xfxQ2(-2, x, mQ2);
        u    = config.pdf -> xfxQ2(1, x, mQ2);
        ubar = config.pdf -> xfxQ2(-1, x, mQ2);
    } else {
        throw std::runtime_error("Unrecognized target type!");
    }

    double s    = config.pdf -> xfxQ2(3, x, mQ2);
    double sbar = config.pdf -> xfxQ2(-3, x, mQ2);
    double b    = config.pdf -> xfxQ2(5, x, mQ2);
    double bbar = config.pdf -> xfxQ2(-5, x, mQ2);

    double c    = config.pdf -> xfxQ2(4, x, mQ2);
    double cbar = config.pdf -> xfxQ2(-4, x, mQ2);
    double t    = config.pdf -> xfxQ2(6, x, mQ2);
    double tbar = config.pdf -> xfxQ2(-6, x, mQ2);

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

double CrossSection::ds_dy_kernel(double k) {
    double x = std::exp(k);

    bool ps_valid = ps.Validate(ENU, x, kernel_y);
    if (!ps_valid) {
        return 0;
    }

    double result = x * ds_dxdy(x, kernel_y);
    return result;
}

double CrossSection::_ds_dy_partonic(double k) {
    double x = std::exp(k);
    return x * ds_dxdy_partonic(x, kernel_y);
}

double CrossSection::ds_dxdy_TMC() {
    return 0.;
}

double CrossSection::ds_dy_TMC() {
    return 0.;
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
    if (config.SF.mass_scheme == "parton") { 
        F.function = &KernelHelper<CrossSection, &CrossSection::_ds_dy_partonic>;
    } else {
        F.function = &KernelHelper<CrossSection, &CrossSection::ds_dy_kernel>;
    }
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

double CrossSection::AlternativeTotalXS(double E) {
    // this is a method using apfelxx/bgr method of integrating (this is their code!)
    double MW2 = config.constants.Mboson2 * SQ(pc->GeV);

    double xl = 1e-9;
    double xu = 1.0 - 1e-9;
    double Ql = sqrt(0.01) * pc->GeV;
    double Qu = 2 * log(500 * MW2);
    double GF2 = pc->GF2;
    double M_target = config.target_mass;

    const double s_tot = SQ(M_target) + 2 * M_target * E;
    // Integration bounds in ln(Q2)
    const double log_Q2min = 2 * log(Ql);
    const double log_Q2max = min(log(s_tot - SQ(M_target)), Qu);
    // const double conv = 0.3894e9;

    const auto dsigmadlnQ2 = [=] (double const& log_Q2) -> double
    {
        const double log_xmin = log_Q2 - log(s_tot - SQ(M_target));
        const double log_xmax = log(xu);

        // Helpers
        const double Q2 = exp(log_Q2);
        // const double Q  = sqrt(Q2);

        const auto dsigmadlnxdlnQ2 = [=] (double const& log_x) -> double
        {
            const double x      = max(exp(log_x), xl);
            const double y      = Q2 / x / ( s_tot - SQ(M_target) );
            const double omy2   = pow(1 - y, 2);
            const double Yplus  = 1 + omy2;
            const double Yminus = 1 - omy2;
            const double W2     = Q2 * (1 - x) / x + SQ(M_target);
            if (W2 < 1.96 * SQ(pc->GeV)) {
                return 0.0;
            }
            double fact = GF2 / 4 / M_PI / x * pow(MW2 / ( MW2 + Q2 ), 2);

            std::array<double, 2> pt{{std::log10(Q2 / SQ(pc->GeV)), std::log10(x)}};

            std::array<int, 2> F1_splc;
            std::array<int, 2> F2_splc;
            std::array<int, 2> F3_splc;

            F1.searchcenters(pt.data(), F1_splc.data());
            F2.searchcenters(pt.data(), F2_splc.data());
            F3.searchcenters(pt.data(), F3_splc.data());

            double f1 = F1.ndsplineeval(pt.data(), F1_splc.data(), 0);
            double f2 = F2.ndsplineeval(pt.data(), F2_splc.data(), 0);
            double xf3 = x * F3.ndsplineeval(pt.data(), F3_splc.data(), 0);
            double fL = f2 - 2*x*f1; // (_F2 * (1.0 + 4.0*SQ(x*m)/Q2) - _FL) / (2. * x);

            return x * fact * ( Yplus * f2
                                - y * y * fL
                                + config.cp_factor * Yminus * xf3 );
        };
        const apfel::Integrator Ixq{dsigmadlnxdlnQ2};
		return Q2 * Ixq.integrate(log_xmin, log_xmax, 1e-5);
    };

    const apfel::Integrator Iq{dsigmadlnQ2};
    const double xsec = Iq.integrate(log_Q2min, log_Q2max, 1e-5);

    return xsec / SQ(pc->cm);
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
    if (config.SF.mass_scheme == "parton") {
        F = { &KernelHelper<CrossSection, &CrossSection::_ds_dxdy_partonic>, dim, this};
    } else {
        F = { &KernelHelper<CrossSection, &CrossSection::ds_dxdy_kernel>, dim, this};
    }
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

// std::tuple<double, double> CrossSection::dsdy_xlims(double s, double y) {
//     double xmin = 0.0;
//     double xmax = 1.0;

//     // Q^2 = s x y > Qmin^2 --> x > Qmin^2 / s y
//     xmin = max(xmin, integral_min_Q2 / (s * y));
//     // W^2 = s y (1 - x) --> x < 1 - Wmin^2 / (2 m E y)
//     xmax = min(xmax, 1.0 - min_W2 / (s * y));
//     // Q^2 = s x y < Qmax^2 --> x < Qmax^2 / s y
//     xmax = min(xmax, integral_max_Q2 / (s * y));

//     return std::make_tuple(xmin, xmax);
// }

// void CrossSection::SetThresholdW2() {

//     double lowest_W2;
//     // If we enable the shallow region, we include the low W^2 region!
//     if (config.XS.enable_shallow_region) {
//         lowest_W2 = SQ(0.938 + 0.13957) * SQ(pc->GeV); // (m_N + m_pi)^2
//     } else {
//         lowest_W2 = 2.0 * SQ(pc->GeV);
//     }

//     min_W2 = lowest_W2;
// }

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
    return calculate_rc_dsdxdy(k, ENU, kernel_x, kernel_y, kernel_L);
}

double CrossSection::calculate_rc_dsdxdy(double z, double E, double x, double y) {
    double M_target = config.target_mass;
    double M_l = config.lepton_mass;
    double s = 2.0 * M_target * E + SQ(M_target);
    double L = log( (s * SQ(1.0 - y + x*y)) / SQ(M_l));
    return calculate_rc_dsdxdy(z, E, x, y, L);
}

double CrossSection::calculate_rc_dsdxdy(double z, double E, double x, double y, double L) {
    double xhat = x * y / (z + y - 1.0);
    double yhat = (z + y - 1.0) / z;
    double zmin = 1.0 - y * (1.0 - x);

    // double M_target = config.target_mass;
    // double M_l = config.lepton_mass;
    // double s = 2.0 * M_target * E + SQ(M_target);
    // double L = log( (s * SQ(1.0 - y + x*y)) / SQ(M_l)); // large logarithm

    std::array<double, 3> pt{{std::log10(E / pc->GeV), std::log10(x), std::log10(y)}};
    std::array<int, 3> xs_splc;
    rc_dsdxdy.searchcenters(pt.data(), xs_splc.data());
    double xs_val = rc_dsdxdy.ndsplineeval(pt.data(), xs_splc.data(), 0);

    double term1 = 0.0;
    if (z > zmin) {
        if (xhat >= 1.0) {
            xhat = 1.0;
        }

        if (yhat >= 1.0) {
            yhat = 1.0;
        }

        if (xhat > 1.0 || yhat > 1.0) {
            std::cout << z << "," << x << "," << y << "," << xhat - 1.0 << "," << yhat - 1.0 << std::endl;   
        }

        std::array<double, 3> pt_hat{{std::log10(E / pc->GeV), std::log10(xhat), std::log10(yhat)}};
        std::array<int, 3> xs_hat_splc;
        rc_dsdxdy.searchcenters(pt_hat.data(), xs_hat_splc.data());
        double xs_hat_val = rc_dsdxdy.ndsplineeval(pt_hat.data(), xs_hat_splc.data(), 0);
        term1 = rc_jacobian(x, y, z) * pow(10.0, xs_hat_val);
    }

    return rc_prefactor * qed_splitting(z) * (term1 - pow(10.0, xs_val));
}

double CrossSection::rc_integrate(double E, double x, double y) {
    Set_Neutrino_Energy(E);
    double M_target = config.target_mass;
    double M_l = config.lepton_mass;
    kernel_y = y;
    kernel_x = x;
    
    double s = 2.0 * M_target * E;// + SQ(M_target);
    kernel_L = log( (s * SQ(1.0 - y + x*y)) / SQ(M_l)); // large logarithm

    double integrate_zmin = 0.0;
    double integrate_zmax = 0.9999;

    if (!ps.Validate(E, x, y)) {
        return 0;
    }

    // gsl_integration_workspace * w = gsl_integration_workspace_alloc(50000);
    gsl_integration_cquad_workspace * w = gsl_integration_cquad_workspace_alloc(2500);
    double result, error;
    size_t neval;

    gsl_function F;
    F.function = &KernelHelper<CrossSection, &CrossSection::rc_kernel>;
    F.params = this;

    // int status = gsl_integration_qag(&F, integrate_zmin, integrate_zmax, 0, 1.e-3, 50000, 6, w, &result, &error);
    int status = gsl_integration_cquad(&F, integrate_zmin, integrate_zmax, 0, 1.e-3, w, &result, &error, &neval);
    if (status != 0) {
        std::cout << "ERR: " << status << std::endl;
    }
    // gsl_integration_workspace_free(w);
    gsl_integration_cquad_workspace_free(w);

    return kernel_L * result;
    
}

void CrossSection::rc_load_dsdxdy(string spline_path) {
    if (rc_dsdxdy_loaded) {
        rc_dsdxdy = photospline::splinetable<>();
    }
    rc_dsdxdy.read_fits(spline_path);
    rc_dsdxdy_loaded = true;
}

}