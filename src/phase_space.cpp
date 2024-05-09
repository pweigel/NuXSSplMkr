#include "NuXSSplMkr/phase_space.h"

namespace nuxssplmkr {

PhaseSpace::PhaseSpace(Configuration& _cfg): config(_cfg) {
    Initialize();
}

void PhaseSpace::Initialize() {
    pc = new nuxssplmkr::PhysConst();
    x_min = config.XS.xmin;
    x_max = config.XS.xmax;
    y_min = 1e-9;
    y_max = 1.0;
    Q2_min = config.XS.Q2min * SQ(pc->GeV); // config
    Q2_max = config.XS.Q2max * SQ(pc->GeV); // config
    W2_min = 4.0  * SQ(pc->GeV); // config/mass
    W2_max = 1e20 * SQ(pc->GeV);

    // If the flag has been set, then change the min W2
    if (config.flag_set_flavor && config.dynamic_W2_min) {
        switch(config.sf_type) {
            case SFType::total:  W2_min = 2.0 * SQ(pc->GeV); break; // TODO
            case SFType::light:  W2_min = 2.0 * SQ(pc->GeV); break; // (m_N + m_pi)^2
            // case SFType::light:  min_W2 = SQ(0.938 + 0.13957) * SQ(pc->GeV); break; // (m_N + m_pi)^2
            // case SFType::charm:  min_W2 = SQ( (0.938 + 1.3) * pc->GeV); break; // (m_N + m_c)^2
            // case SFType::bottom: min_W2 = SQ( (0.938 + 4.5) * pc->GeV); break; // (m_N + m_b)^2
            // case SFType::top:    min_W2 = SQ( (0.938 + 173.0) * pc->GeV); break; // (m_N + m_t)^2, TODO: get right val
            case SFType::charm:  W2_min = SQ( (0.938 + 1.870) * pc->GeV); break; // (m_N + m_D)^2
            case SFType::bottom: W2_min = SQ( (0.938 + 5.279) * pc->GeV); break; // (m_N + m_B)^2
            case SFType::top:    W2_min = SQ( (0.938 + 173.0) * pc->GeV); break; // (m_N + m_t)^2, TODO: get right val
            default:             W2_min = 2.0 * SQ(pc->GeV); break; // TODO
        }
    } else if (config.dynamic_W2_min) {
        std::cout << "WARNING: config flavor type has not been set when initializing the phase space and dynamic W2min was enabled!!" << std::endl;
    }
}

void PhaseSpace::SetFinalStateMass(double m_fs) {
    /* Here, we compute the constraints on the phase space given
    the conditions specified in the configuration and final state
    constraints */

    // Base constraints
    // xmin = config[]

    // Final state
    // W2_min = pow(pc->m_N + m_fs, 2);
}

bool PhaseSpace::Validate(double E) {
    /*
    TODO
    */
    flag=0;
    double s = 2 * pc->m_N * E + pow(pc->m_N, 2);

    if (s < Q2_min || s < W2_min) {
        flag = 1;
        return false;
    }

    return true;
}

bool PhaseSpace::Validate(double E, double y) {
    /*
    TODO
    */
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

    double s = 2 * pc->m_N * E + pow(pc->m_N, 2);
    double Q2 = (s - pow(pc->m_N, 2)) * x * y;
    double W2 = s * y * (1.0 - x);

    if (debug) {
        std::cout << "s = " << s / SQ(pc->GeV) << std::endl;
        std::cout << "Q2 = " << Q2 / SQ(pc->GeV) << std::endl;
        std::cout << "W2 = " << W2 / SQ(pc->GeV) << std::endl;
    }
    // double eta = 1.0; // Q^2 / m^2

    if ((Q2 < Q2_min) || (Q2 > Q2_max) || (s > Q2_max)) {
        /* 
        Constraints from:
        PDF minimum/maximum Q^2
        If CKMT/PCAC is on, it is lowered to some value like ~0.1 GeV^2
        */
        flag = 1;
        return false;
    }

    if ((x < x_min) || (x > x_max)) {
        /* 
        Constraints from:
        PDF minimum x

        */
        flag = 2;
        return false;
    }

    if ((y < y_min) || (y > y_max)) {
        /*
        Constraints from:

        */
        flag = 3;
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
    std::cout << "Phase space constraints: " << std::endl;
    std::cout << "  (xmin, xmax) = " << x_min << ", " << x_max << std::endl;
    std::cout << "  (ymin, ymax) = " << y_min << ", " << y_max << std::endl;
    std::cout << "  (Q2min, Q2max) [GeV^2] = " << Q2_min / SQ(pc->GeV) << ", " << Q2_max / SQ(pc->GeV)<< std::endl;
    std::cout << "  (W2min, W2max) [GeV^2] = " << W2_min / SQ(pc->GeV) << ", " << W2_max / SQ(pc->GeV) << std::endl;
}
// bool CrossSection::PhaseSpaceIsGood_Q2(double x, double Q2, double E) {
//     Set_Neutrino_Energy(E);
//     double s = 2.0 * M_iso * E + SQ(M_iso);

//     // First check that the x is within the bounds of the SF grids  
//     if ((x < config.SF.xmin) || (x > config.SF.xmax)) {
//         return false;
//     }

//     // Check that the Q2 is within the bounds of the SF grids
//     if ( (Q2 < (config.SF.Q2min * SQ(pc->GeV))) || (Q2 > (config.SF.Q2max * SQ(pc->GeV))) ) {
//         return false;
//     }
    
//     // integral_max_Q2 = s;
//     // Check that the Q2 is within the integration bounds
//     if ( ( Q2 < integral_min_Q2 ) || (Q2 > integral_max_Q2) ) {
//         return false;
//     }

//     // Make sure y is less than 1
//     double y = Q2 / (s * x);
//     if (y > 1) {
//         return false;
//     }
    
//     // Calculate W^2 = Q^2 (1/x - 1) + M_N^2
//     double W2 =  Q2 * (1.0 - x) / (x + 1e-15) + SQ(M_iso); // TODO: target mass

//     // try to fix top? this is the fake threshold from the pdfs
//     if ( (config.sf_type == SFType::top) && ( Q2 < 34 * SQ(pc->GeV))) {
//         return false;
//     }

//     // Check W^2 threshold
//     if ( W2 < min_W2 ) {
//         return false;
//     }

//     return true;
// }

// bool CrossSection::PhaseSpaceIsGood(double x, double y, double E) {
//     double s = 2.0 * M_iso * E + SQ(M_iso);
//     double Q2 = (s - SQ(M_iso)) * x * y; 
//     return PhaseSpaceIsGood_Q2(x, Q2, E);
// }

}