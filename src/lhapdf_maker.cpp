#include "NuXSSplMkr/lhapdf_maker.h"

namespace nuxssplmkr {

LHAPDFMaker::LHAPDFMaker(Configuration &_config)
    : config(_config)
{ 
    suffix = "";
    switch(config.mode) {
        case 1: suffix = ""; break;
        case 2: suffix = "_TMC"; break;
        case 3: suffix = "_CKMT"; break;
        case 4: suffix = "_PCAC"; break;
    }
    const char* lha_env_var = std::getenv("LHAPDF_DATA_PATH")    ;
    const std::string lha_path = lha_env_var == NULL ? std::string(".") : std::string(lha_env_var);
    pdf_path = lha_path + "/" + config.general.unique_name + suffix +"_SF/";

    std::cout << "Saving structure functions to: " << std::endl;
    std::cout << "  " << pdf_path << std::endl;

    boost::filesystem::path out_folder = pdf_path;
    if (!boost::filesystem::exists(out_folder)) {
        boost::filesystem::create_directories(out_folder);
    }
}

std::vector<std::string> LHAPDFMaker::MakeSet(string datapath) {
    std::string sfs[3] = {"F2", "F1", "F3"};
    std::string flavors[4] = {"light", "charm", "bottom", "top"};
    std::string projectiles[2] = {"neutrino", "antineutrino"};

    std::unordered_map<std::string, std::vector<std::vector<double>>> sf_data;
    std::vector<std::string> sf_codes;

    bool initialized_xQ2 = false;
    std::vector<double> logx_values;
    std::vector<double> logQ2_values;
    unsigned int nQ2 = 1;
    unsigned int nx = 1;
    for (auto &_proj : projectiles) {
        for (auto &_sf : sfs) {
            for (auto &_fl : flavors) {
                std::string code = SF_INTERACTION_CODE["CC"] + SF_PARTICLE_CODE[_proj] + SF_NUMBER_CODES[_sf] + SF_FLAVOR_CODES[_fl];

                std::string infile = datapath + "/" + _sf + "_" + _proj + "_" + config.target + "_" + _fl + suffix + ".grid";
                std::cout << infile << std::endl;
                std::ifstream gridfile(infile);
                std::string line;

                if (!initialized_xQ2) {
                    // The first line contains the log Q2 values
                    std::getline(gridfile, line);
                    std::istringstream logQ2_stream(line);
                    double _logQ2;
                    while (logQ2_stream >> _logQ2) {logQ2_values.push_back(_logQ2);}

                    // The second line contains the log x values
                    std::getline(gridfile, line);
                    std::istringstream logx_stream(line);
                    double _logx;
                    while (logx_stream >> _logx) {logx_values.push_back(_logx);}

                    nQ2 = logQ2_values.size();
                    nx = logx_values.size();

                    initialized_xQ2 = true;
                } else {
                    // 2 lines of header we need to skip
                    std::getline(gridfile, line);
                    std::getline(gridfile, line);
                }

                // create data
                std::vector<std::vector<double>> data(nQ2, std::vector<double>(nx, 0.0));

                for (unsigned int i = 0; i < nQ2; i++){
                    std::getline(gridfile, line);
                    std::istringstream linestream(line);
                    std::string val;
                    unsigned int j = 0;
                    while (std::getline(linestream, val, ',')) {
                        data[i][j] =  std::stod(val);
                        j += 1;
                    }
                }
                gridfile.close();

                sf_data[code] = data;
                sf_codes.push_back(code);
            }
        }
    }

    std::ofstream outfile(pdf_path + config.general.unique_name+suffix+"_SF_0000.dat");

    outfile << "PdfType: central\nFormat: lhagrid1\n---\n";
    outfile << std::scientific << std::setprecision(11) << "   ";
    for (unsigned int i = 0; i < nx; i++) {
        double x = std::pow(10, logx_values.at(i));
        outfile << x << "   ";
    }
    outfile << std::endl << "   ";
    for (unsigned int i = 0; i < nQ2; i++) {
        double Q2 = std::pow(10, logQ2_values.at(i));
        outfile << std::sqrt(Q2) << "   ";
    }
    outfile << std::endl << "   " << std::setw(2);
    for (auto &sf_code : sf_codes) {
        outfile << sf_code << "   ";
    }
    outfile << std::endl;

    outfile << std::scientific << std::setprecision(8);
    for (unsigned int j = 0; j < nx; j++) {
        for (unsigned int i = 0; i < nQ2; i++) {
            outfile << "   "; // << logQ2_values.at(i) << " " << logx_values.at(j) << " ";
            for (auto &sf_code : sf_codes) {
                outfile << sf_data[sf_code][i][j] << "   ";
            }
            outfile << std::endl;
        }
    }
    outfile << "---";
    outfile.close();

    return sf_codes;
}

void LHAPDFMaker::MakeInfo(std::vector<std::string> codes) {
    // The following code is from the nnpdf collab code!
    // LHAPDF6 HEADER
    std::ofstream infodata;
    infodata.open(pdf_path+config.general.unique_name + suffix + "_SF"+".info");
    int nrep = 1; // TODO!
    infodata << "SetDesc: \"Structure functions!\"" << endl;
    infodata << "SetIndex: " << endl;
    infodata << "Authors: Philip Weigel" << endl;
    infodata << "Reference: N/A" << endl;
    infodata << "Format: lhagrid1" << endl;
    infodata << "DataVersion: 1" << endl;
    infodata << "NumMembers: " << nrep << endl;
    infodata << "Particle: 2212" << endl;
    infodata << "Flavors: [";
    for (int i = 0; i < codes.size()-1; i++) {
        infodata << codes.at(i) << ", ";
    }
    infodata << codes.at(codes.size()-1) << "]" << std::endl;
    infodata << "OrderQCD: " << config.SF.pto << endl;
    infodata << "FlavorScheme: variable" << endl; // TODO: currently always set as variable which isn't always the case
    infodata << "NumFlavors: 6" << endl; // should be a part of the config
    infodata << "ErrorType: replicas" << endl; // get from lhapdf

    infodata.precision(7);
    infodata << scientific;
    infodata << "XMin: "<< config.SF.xmin << endl;
    infodata << "XMax: "<< config.SF.xmax << endl;
    infodata << "QMin: "<< to_string(sqrt(config.SF.Q2min)) << endl;
    infodata << "QMax: "<< to_string(sqrt(config.SF.Q2max)) << endl;
    infodata << "MZ: "  << config.constants.MassZ << endl;
    infodata << "MUp: 0\nMDown: 0\nMStrange: 0" << std::endl;
    infodata << "MCharm: "  << config.pdf_info.pdf_quark_masses[4] << endl;
    infodata << "MBottom: " << config.pdf_info.pdf_quark_masses[5] << endl;
    infodata << "MTop: "    << config.pdf_info.pdf_quark_masses[6] << endl;

    LHAPDF::PDFInfo info = config.pdf->info();
    if (info.has_key("AlphaS_MZ")) {
        std::string AlphaS_MZ = info.get_entry_as<string>("AlphaS_MZ");
        infodata << "AlphaS_MZ: " << AlphaS_MZ << endl;
    }
    if (info.has_key("AlphaS_OrderQCD")) {
        std::string AlphaS_OrderQCD = info.get_entry_as<string>("AlphaS_OrderQCD");
        infodata << "AlphaS_OrderQCD: " << AlphaS_OrderQCD << endl;
    }
    if (info.has_key("AlphaS_Type")) {
        std::string AlphaS_Type = info.get_entry_as<string>("AlphaS_Type");
        infodata << "AlphaS_Type: " << AlphaS_Type << endl;
    }

    if (info.has_key("AlphaS_Qs")) {
        std::vector<double> AlphaS_Qs = info.get_entry_as<vector<double>>("AlphaS_Qs");
        infodata << "AlphaS_Qs: [";
        for (int i = 0; i < AlphaS_Qs.size()-1; i++) {
            infodata << AlphaS_Qs.at(i) << ", ";
        }
        infodata << AlphaS_Qs.at(AlphaS_Qs.size()-1) << "]" << endl;
    }
    if (info.has_key("AlphaS_Vals")) {
        std::vector<double> AlphaS_Vals = info.get_entry_as<vector<double>>("AlphaS_Vals");
        infodata << "AlphaS_Vals: [";
        for (int i = 0; i < AlphaS_Vals.size()-1; i++) {
            infodata << AlphaS_Vals.at(i) << ", ";
        }
        infodata << AlphaS_Vals.at(AlphaS_Vals.size()-1) << "]" << endl;
    }
    // infodata << "AlphaS_Lambda4: 0.342207" << endl;
    // infodata << "AlphaS_Lambda5: 0.239" << endl;
    infodata.close();
}

}
