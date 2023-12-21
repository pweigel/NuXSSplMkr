#include "cross_section.h"

namespace nuxssplmkr {

CrossSection::CrossSection(StructureFunction& _sf): sf(_sf) {
  
}

double CrossSection::ds_dxdy(double x, double y, double Q2) {
    double _FL = sf.FL(x, Q2);
    double _F2 = sf.F2(x, Q2);
    double _F1 = (_F2 - _FL) / (2. * x);
    double _F3 = sf.F3(x, Q2) / x;
    return 0.;
    // return x*y*y * _F1 + (1 - y) * _F2 + sf.CP_Factor * x * y * (1 - y/2) * _F3;
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

}