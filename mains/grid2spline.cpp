#include "NuXSSplMkr/configuration.h"
#include "NuXSSplMkr/physconst.h"
#include "NuXSSplMkr/structure_function.h"
#include "NuXSSplMkr/spline_maker.h"
#include <boost/filesystem.hpp>

int main(int argc, char* argv[]){
    const string infile = argv[1];
    const string outfile = argv[2];

    nuxssplmkr::SplineMaker splmkr = nuxssplmkr::SplineMaker();
    splmkr.MakeSpline(infile, outfile);

}