#include "NuXSSplMkr/phase_space.h"

namespace nuxssplmkr {

PhaseSpace::PhaseSpace(Configuration& _cfg): config(_cfg) {
    Initialize();
}

void PhaseSpace::Initialize() {
    pc = new nuxssplmkr::PhysConst();
    x_min = config.XS.xmin;
    x_max = config.XS.xmax;
    y_min = config.XS.ymin;
    y_max = config.XS.ymax;
    Q2_min = config.XS.Q2min * SQ(pc->GeV);
    Q2_max = config.XS.Q2max * SQ(pc->GeV);
    W2_min = config.XS.W2min * SQ(pc->GeV);
    W2_max = 1e20 * SQ(pc->GeV);

    debug = config.general.debug;
}

void PhaseSpace::SetFinalStateMass(double m_fs) {
    /* Here, we compute the constraints on the phase space given
    the conditions specified in the configuration and final state
    constraints */
}

bool PhaseSpace::Validate(double E) {
    flag=0;
    double s = 2 * pc->m_N * E + pow(pc->m_N, 2);

    if (s < Q2_min || s < W2_min) {
        flag = 1;
        return false;
    }

    return true;
}

bool PhaseSpace::Validate(double E, double y) {
    flag = 0;
    Validate(E);
    
    return true;
}

bool PhaseSpace::Validate_xQ2(double E, double x, double Q2) {
    // E and Q^2 are given in units of eV and eV^2
    // std::cout << E << ", " << x << ", " << Q2 << std::endl;
    if ((Q2 < Q2_min) || (Q2 > 2.0 * pc->m_N * E)) {
        return false;
    }
    double _xmin = max(x_min, Q2 / (2*pc->m_N*E));
    if ((x < _xmin) || (x >= 1)) {
        return false;
    }
    double W2 = Q2 * (1 - x) / x + pow(pc->m_N, 2);
    double _W2min = 4.0 * pow(pc->GeV, 2);
    if (W2 < _W2min) {
        return false;
    }
    return true;
}

bool PhaseSpace::Validate(double E, double x, double y) {
    /* 
    Check that the point (E, x, y) is valid given the constraints
    on W2 and Q2. This should be called at every point during the
    integrations.
    */
    flag = 0; // No problem!

    double s = 2 * pc->m_N * E + pow(pc->m_N, 2);
    double Q2 = (s - pow(pc->m_N, 2)) * x * y;
    // double W2 = s * y * (1.0 - x);
    double W2 = Q2 * (1.0 - x) / x + pow(pc->m_N, 2);

    // use other ps check
    double M_l = 0.1057 * pc->GeV;
    double _xmin = max(x_min, SQ(M_l) / (2.0 * pc->m_N * (E - M_l)));
    // double _xmin = max(x_min, Q2 / (2.0 * pc->m_N * E));
    // double _xmin = max(x_min, Q2_min / (s * y));
    double _xmax = min(x_max, 1.0);

    double a = (1.0 - SQ(M_l) * (1.0 / (2.0 * pc->m_N * E * x) + 1.0 / (2 * SQ(E)))) / (2.0 + pc->m_N * x / E);
    double b = sqrt(SQ(1.0 - SQ(M_l) / (2.0 * pc->m_N * E * x)) - SQ(M_l / E)) / (2.0 + pc->m_N * x / E);
    double _ymin = a - b;
    double _ymax = a + b;

    if (debug) {
        std::cout << "_xmin = " << _xmin << std::endl;
        std::cout << "_xmax = " << _xmax << std::endl;
        std::cout << "_ymin = " << _ymin << std::endl;
        std::cout << "_ymax = " << _ymax << std::endl;
    }

    if ((x < _xmin) || (x > _xmax) || (y < _ymin) || (y > _ymax)) {
        flag = 1;
        return false;
    }

    if (debug) {
        std::cout << "s = " << s << std::endl;
        std::cout << "Q2 = " << Q2 << std::endl;
        std::cout << "W2 = " << W2 << std::endl;
    }
    // double eta = 1.0; // Q^2 / m^2

    if ((Q2 < Q2_min) || (Q2 > s) || (Q2_max > Q2_max)) {
        /* 
        Constraints from:
        PDF minimum/maximum Q^2
        If CKMT/PCAC is on, it is lowered to some value like ~0.1 GeV^2
        */
        flag = 1;
        return false;
    }

    // if ((x < x_min) || (x > x_max)) {
    //     /* 
    //     Constraints from:
    //     PDF minimum x

    //     */
    //     flag = 2;
    //     return false;
    // }

    // if ((y < y_min) || (y > y_max)) {
    //     /*
    //     Constraints from:

    //     */
    //     flag = 3;
    //     return false;
    // }

    if ((W2 < W2_min) || (W2 > W2_max)) {
        /*
        Constraints from:
        Definition of "deep": W^2 > 4 (sometimes 3.5?)
        */
        flag = 4;
        return false;
    }

    return true;
}


void PhaseSpace::Print() {
    std::cout << "Phase space limits: " << std::endl;
    std::cout << "  (xmin, xmax) = " << x_min << ", " << x_max << std::endl;
    std::cout << "  (ymin, ymax) = " << y_min << ", " << y_max << std::endl;
    std::cout << "  (Q2min, Q2max) [GeV^2] = " << Q2_min / SQ(pc->GeV) << ", " << Q2_max / SQ(pc->GeV)<< std::endl;
    std::cout << "  (W2min, W2max) [GeV^2] = " << W2_min / SQ(pc->GeV) << ", " << W2_max / SQ(pc->GeV) << std::endl;
}

}