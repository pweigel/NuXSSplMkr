#include "phase_space.h"

PhaseSpace::PhaseSpace() {

}

PhaseSpace::Compute() {
    /* Here, we compute the constraints on the phase space given
    the conditions specified in the configuration and final state
    constraints */

    // Final state

}

PhaseSpace::Validate(double E, double x, double y) {
    /* 
    Check that the point (E, x, y) is valid given the constraints
    on W2 and Q2. This should be called at every point during the
    integrations.
    */
    double s = 2 * m_N * E + pow(m_N, 2);
    double Q2 = (s - pow(m_N, 2)) * x * y;
    double W2 = s * y * (1.0 - x);
    double eta = 1.0; // Q^2 / m^2

    if ((Q2 < Q2_min) || (Q2 > Q2max)) {
        /* 
        Constraints from:
        PDF minimum/maximum Q^2
        If CKMT/PCAC is on, it is lowered to some value like ~0.1 GeV^2
        */
        return false;
    }

    if ((x < x_min) || (x > x_max)) {
        /* 
        Constraints from:
        PDF minimum x

        */
        return false;
    }

    if ((y < y_min) || (y > y_max)) {
        /*
        Constraints from:

        */
        return false;
    }

    if ((W2 < W2_min) || (W2 > W2_max)) {
        /*
        Constraints from:
        Definition of "deep": W^2 > 4 (sometimes 3.5?)
        */
        return false;
    }

    return true;
}

bool CrossSection::PhaseSpaceIsGood_Q2(double x, double Q2, double E) {
    Set_Neutrino_Energy(E);
    double s = 2.0 * M_iso * E + SQ(M_iso);

    // First check that the x is within the bounds of the SF grids  
    if ((x < config.SF.xmin) || (x > config.SF.xmax)) {
        return false;
    }

    // Check that the Q2 is within the bounds of the SF grids
    if ( (Q2 < (config.SF.Q2min * SQ(pc->GeV))) || (Q2 > (config.SF.Q2max * SQ(pc->GeV))) ) {
        return false;
    }
    
    // integral_max_Q2 = s;
    // Check that the Q2 is within the integration bounds
    if ( ( Q2 < integral_min_Q2 ) || (Q2 > integral_max_Q2) ) {
        return false;
    }

    // Make sure y is less than 1
    double y = Q2 / (s * x);
    if (y > 1) {
        return false;
    }
    
    // Calculate W^2 = Q^2 (1/x - 1) + M_N^2
    double W2 =  Q2 * (1.0 - x) / (x + 1e-15) + SQ(M_iso); // TODO: target mass

    // try to fix top? this is the fake threshold from the pdfs
    if ( (config.sf_type == SFType::top) && ( Q2 < 34 * SQ(pc->GeV))) {
        return false;
    }

    // Check W^2 threshold
    if ( W2 < min_W2 ) {
        return false;
    }

    return true;
}

bool CrossSection::PhaseSpaceIsGood(double x, double y, double E) {
    double s = 2.0 * M_iso * E + SQ(M_iso);
    double Q2 = (s - SQ(M_iso)) * x * y; 
    return PhaseSpaceIsGood_Q2(x, Q2, E);
}