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
    if ((config.mode == 3) || (config.mode == 4)) { // we do this in the config now
        Q2_min = 0.01 * SQ(pc->GeV);
    }
    Q2_max = config.XS.Q2max * SQ(pc->GeV);
    W2_min = config.XS.W2min * SQ(pc->GeV);
    if ((config.mode == 3) || (config.mode == 4)) {
        W2_min = SQ(0.938 + 0.140) * SQ(pc->GeV);
    }
    W2_max = 1e20 * SQ(pc->GeV);

    debug = config.general.debug;
}

void PhaseSpace::SetFinalStateMass(double m_fs) {
    /* Here, we compute the constraints on the phase space given
    the conditions specified in the configuration and final state
    constraints */
}

bool PhaseSpace::Validate(double E) {
    double M_target = config.target_mass;
    double s = 2 * M_target * E + SQ(M_target);

    flag = 0;
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


bool PhaseSpace::Validate(double E, double x, double y) {
    /* 
    Check that the point (E, x, y) is valid given the constraints
    on W2 and Q2. This should be called at every point during the
    integrations.
    */
    flag = 0; // No problem!
    double M_target = config.target_mass;
    double M_l = config.lepton_mass;

    double s = 2 * M_target * E + pow(M_target, 2);
    double Q2 = (s - pow(M_target, 2)) * x * y;
    double W2 = Q2 * (1.0 - x) / x + pow(M_target, 2);

    // use other ps check
    double _xmin = max(x_min, SQ(M_l) / (2.0 * M_target * (E - M_l)));
    double _xmax = min(x_max, 1.0);

    double a = (1.0 - SQ(M_l) * (1.0 / (2.0 * M_target * E * x) + 1.0 / (2 * SQ(E)))) / (2.0 + M_target * x / E);
    double b = sqrt(SQ(1.0 - SQ(M_l) / (2.0 * M_target * E * x)) - SQ(M_l / E)) / (2.0 + M_target * x / E);
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

    if ((Q2 < Q2_min) || (Q2 > s - SQ(M_target)) || (Q2_max > Q2_max)) {
        /* 
        Constraints from:
        PDF minimum/maximum Q^2
        If CKMT/PCAC is on, it is lowered to some value like ~0.1 GeV^2
        */
        flag = 1;
        return false;
    }

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
    std::cout << "Phase space limits from config: " << std::endl;
    std::cout << "  (xmin, xmax) = " << x_min << ", " << x_max << std::endl;
    std::cout << "  (ymin, ymax) = " << y_min << ", " << y_max << std::endl;
    std::cout << "  (Q2min, Q2max) [GeV^2] = " << Q2_min / SQ(pc->GeV) << ", " << Q2_max / SQ(pc->GeV)<< std::endl;
    std::cout << "  (W2min, W2max) [GeV^2] = " << W2_min / SQ(pc->GeV) << ", " << W2_max / SQ(pc->GeV) << std::endl;
}

}