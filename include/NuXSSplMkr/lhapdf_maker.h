#ifndef __LHAPDF_MAKER_H
#define __LHAPDF_MAKER_H

#include "NuXSSplMkr/configuration.h"
#include "NuXSSplMkr/tools.h"
#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/GridPDF.h"

namespace nuxssplmkr {
  
static unordered_map<string, string> SF_INTERACTION_CODE {
    {"CC", "1"},
    {"NC", "2"},
};

static unordered_map<string, string> SF_PARTICLE_CODE {
    {"neutrino", "0"},
    {"antineutrino", "1"},
};

static unordered_map<string, string> SF_NUMBER_CODES {
    {"F2", "1"},
    {"FL", "2"},
    {"xF3", "3'"},
    {"F1", "4"},
    {"F3", "5"},
};

static unordered_map<string, string> SF_FLAVOR_CODES {
    {"total", "0"},
    {"light", "1"},
    {"charm", "2"},
    {"bottom", "3"},
    {"top", "4"},
};

class LHAPDFMaker {
    private:
        Configuration config;
    public:
        LHAPDFMaker(Configuration& _config);
        ~LHAPDFMaker() {};

        std::vector<std::string> MakeSet(string datapath);
        void MakeInfo(std::vector<std::string> codes);
};

}

#endif