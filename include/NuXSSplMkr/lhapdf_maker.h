#ifndef __LHAPDF_MAKER_H
#define __LHAPDF_MAKER_H

#include "NuXSSplMkr/configuration.h"
#include "NuXSSplMkr/tools.h"
#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/GridPDF.h"
#include <boost/filesystem.hpp>

namespace nuxssplmkr {

class LHAPDFMaker {
    private:
        Configuration config;
        string suffix;
        string pdf_path;
    public:
        LHAPDFMaker(Configuration& _config);
        ~LHAPDFMaker() {};

        std::vector<std::string> MakeSet(string datapath);
        void MakeInfo(std::vector<std::string> codes);
};

}

#endif