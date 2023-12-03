#include "structure_function.h"

namespace nuxssplmkr {

StructureFunction::StructureFunction(Configuration &_config)
  : config(_config)
{
    pc = new PhysConst();
}

void InitializeAPFEL() {
    
}

double StructureFunction::F1(double x, double Q2) {
    return 0.;
}

double StructureFunction::F2(double x, double Q2) {
    return 0.;
}

double StructureFunction::xF3(double x, double Q2) {
    return 0.;
}

double StructureFunction::F3(double x, double Q2) {
    return 0.;
}

}