#ifndef __SPLINE_MAKER_H
#define __SPLINE_MAKER_H

#include "NuXSSplMkr/configuration.h"
#include "NuXSSplMkr/tools.h"
#include <photospline/splinetable.h>
#include <deque>

namespace nuxssplmkr {
  
class SplineMaker {
    private:

    public:
        SplineMaker();
        ~SplineMaker() {};

        void MakeSpline(string infile, string outfile);
};

}

#endif