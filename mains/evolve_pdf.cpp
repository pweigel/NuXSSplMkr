#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cmath>
#include <boost/filesystem.hpp>
#include "APFEL/APFEL.h"
#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/GridPDF.h"

using namespace std;

int main(int argc, char* argv[]){
    if (argc != 5) {
        std::cout << "Not enough/too many inputs!" << std::endl;
        std::cout << "Usage: evolve_pdfs INPDF OUTPDF REPLICA ENABLE_SMALLX" << std::endl;
        return 1;
    }

    std::string inpdf = argv[1];  // input pdfset name
    std::string outpdf = argv[2]; // output pdfset name
    std::string replica = argv[3]; // replica number
    bool enable_small_x = !!(std::stoi(argv[4])); // use small-x resummation?
    int rep = std::stoi(replica);

    auto replica_string = std::string(4 - std::min(4, static_cast<int>(replica.length()) ), '0') + replica;
    boost::filesystem::path out_folder = outpdf+"/";
    if (!boost::filesystem::exists(out_folder)) {
        boost::filesystem::create_directories(out_folder);
    }

    LHAPDF::PDF* pdf = LHAPDF::mkPDF(inpdf, rep);
    LHAPDF::PDFInfo info = pdf->info();

    std::ofstream outfile;
    outfile.open(outpdf + "/" + outpdf + "_" + replica_string + ".dat");

    if (replica == "0") outfile << "PdfType: central" << std::endl;
    else outfile << "PdfType: replica" << std::endl;
    outfile << "Format: lhagrid1" << std::endl << "---" << std::endl;

    int pto = info.get_entry_as<int>("OrderQCD");
    double xmin = info.get_entry_as<double>("XMin");
    double Q0 = info.get_entry_as<double>("QMin");
    if (enable_small_x) {
        std::cout << "Small-x resummation is enabled, setting Q0 = sqrt(3.5) GeV!" << std::endl;
        Q0 = std::sqrt(3.5);
    }
    double Q_max = 1e5;

    double MZ, MUp, MDown, MStrange, MCharm, MBottom, MTop, AlphaS_MZ;
    if (info.has_key("MZ")) MZ = info.get_entry_as<double>("MZ");
        else MZ = 91.1870;
    if (info.has_key("MUp")) MUp = info.get_entry_as<double>("MUp");
        else MUp = 0.0;
    if (info.has_key("MDown")) MDown = info.get_entry_as<double>("MDown");
        else MDown = 0.0;
    if (info.has_key("MStrange")) MStrange = info.get_entry_as<double>("MStrange");
        else MStrange = 0.0;
    if (info.has_key("MCharm")) MCharm = info.get_entry_as<double>("MCharm");
        else MCharm = 1.3;
    if (info.has_key("MBottom")) MBottom = info.get_entry_as<double>("MBottom");
        else MBottom = 4.75;
    if (info.has_key("MTop")) MTop = info.get_entry_as<double>("MTop");
        else MTop = 172.0;
    if (info.has_key("AlphaS_MZ")) AlphaS_MZ = info.get_entry_as<double>("AlphaS_MZ");
        else AlphaS_MZ = 0.118;

    APFEL::SetPDFSet(inpdf);
    APFEL::SetReplica(rep);
    APFEL::SetTheory("QCD");
    APFEL::SetQLimits(Q0, Q_max);
    APFEL::SetPerturbativeOrder(pto);
    APFEL::SetPDFEvolution("exactalpha");
    APFEL::SetFastEvolution(false);
    APFEL::SetAlphaEvolution("exact");
    APFEL::SetVFNS();
    APFEL::SetMaxFlavourPDFs(6);
    APFEL::SetMaxFlavourAlpha(6);
    APFEL::SetNumberOfGrids(3);
    APFEL::SetZMass(MZ);
    APFEL::SetGridParameters(1, 130, 3, xmin);
    APFEL::SetGridParameters(2, 60, 5, 1e-1);
    APFEL::SetGridParameters(3, 20, 5, 8e-1);
    APFEL::SetPoleMasses(MCharm, MBottom, MTop);
    APFEL::SetAlphaQCDRef(AlphaS_MZ, MZ);

    int nQ2 = 100;
    int nlogx = 100;
    int nlinx = 50;
    int nx = nlogx + nlinx;

    double logQ2_min = std::log10(Q0*Q0);
    double logQ2_max = std::log10(Q_max*Q_max);
    double logx_min = std::log10(xmin);
    double logx_max = -1;
    double x_max = 1;

    double Q2_array[nQ2];
    double x_array[nx];
    double dlogx = (logx_max - logx_min) / nlogx;
    double dx = (x_max - pow(10, logx_max)) / (nlinx - 1);
    outfile << std::scientific << std::setprecision(11) << "  ";
    for (int i = 0; i < nlogx; i++) {
        double x = pow(10, logx_min + dlogx * i);
        x_array[i] = x;
        outfile << x << " ";
    }
    for (int i = 0; i < nlinx; i++) {
        double x = pow(10, logx_max) + dx * i;
        x_array[nlogx+i] = x;
        outfile << x << " ";
    }
    outfile << std::endl << "  ";

    double dlogQ2 = (logQ2_max - logQ2_min) / (nQ2 - 1);
    for (int i = 0; i < nQ2; i++) {
        double Q2 = pow(10, logQ2_min + dlogQ2 * i);
        Q2_array[i] = Q2;
        outfile << sqrt(Q2) << " ";
    }
    outfile << std::endl << "  " << std::setw(2);
    for (int p = -6; p < 7; p++) {
        if (p == 0) {
            outfile << "21 ";
        } else{
            outfile << p << " ";
        }
    }
    outfile << std::endl;

    APFEL::SetSmallxResummation(enable_small_x, "NLL");
    APFEL::InitializeAPFEL();

    double evolved_pdf[nQ2][nx][13];
    for (int i = 0; i < nQ2; i++) {
        APFEL::EvolveAPFEL(Q0, std::sqrt(Q2_array[i]));
        // std::cout << Q2_array[i] << std::endl;

        for (int j = 0; j < nx; j++) {
            double x = x_array[j];
            for (int p = -6; p < 7; p++) {
                double xf = APFEL::xPDF(p, x);
                evolved_pdf[i][j][p+6] = xf;
            }
        }
    }

    outfile << std::scientific << std::setprecision(8);
    for (int j = 0; j < nx; j++) {
        for (int i = 0; i < nQ2; i++) {
            // double x = x_array[j];
            // std::cout << std::scientific << std::setprecision(1) << x << " " << std::setprecision(4);
            outfile << " ";
            for (int p = -6; p < 7; p++) {
                outfile << evolved_pdf[i][j][p+6] << " ";
            }
            outfile << std::endl;
        }
    }
    outfile << "---";
    outfile.close();

    if (!boost::filesystem::exists(outpdf + "/" + outpdf + ".info")) {
        // Get some useful things
        int nreps = info.get_entry_as<int>("NumMembers");
        std::string error_type = info.get_entry_as<std::string>("ErrorType");
        double error_conf_level = info.get_entry_as<double>("ErrorConfLevel");

        // The following code is a modified version of the nnpdf collab code!
        // LHAPDF6 HEADER
        std::ofstream infodata;
        infodata.open(outpdf + "/" + outpdf + ".info");
        infodata << "SetDesc: \"Modified version of " << inpdf <<  ", not for public use!\"" << endl;
        infodata << "SetIndex: " << endl;
        infodata << "Authors: Philip Weigel" << endl;
        infodata << "Reference: N/A" << endl;
        infodata << "Format: lhagrid1" << endl;
        infodata << "DataVersion: 1" << endl;
        infodata << "NumMembers: " << nreps << endl;
        infodata << "Flavors: [";
        for (int p = -6; p < 7; p++) {
            if      (p == 0) infodata << "21, ";
            else if (p == 6) infodata << "6";
            else             infodata << p << ", ";
        }
        infodata << "]" << endl;
        infodata << "OrderQCD: " << APFEL::GetPerturbativeOrder() << endl;
        infodata << "FlavorScheme: variable" << endl;
        infodata << "NumFlavors: " << APFEL::GetMaxFlavourPDFs() << endl;
        infodata << "ErrorType: " << error_type << endl;
        infodata << "ErrorConfLevel: " << error_conf_level << endl;

        infodata.precision(7);
        infodata << scientific;
        infodata << "XMin: "<< x_array[0] << endl;
        infodata << "XMax: "<< x_array[nx-1] << endl;
        infodata << "QMin: "<< Q0 << endl;
        infodata << "QMax: "<< Q_max << endl;
        infodata << "MZ: "  << APFEL::GetZMass() << endl;
        infodata << "MUp: " << MUp << endl;
        infodata << "MDown: " << MDown << endl;
        infodata << "MStrange: " << MStrange << endl;
        infodata << "MCharm: "  << APFEL::GetThreshold(4) << endl;
        infodata << "MBottom: " << APFEL::GetThreshold(5) << endl;
        infodata << "MTop: "    << APFEL::GetThreshold(6) << endl;
        infodata << fixed << "AlphaS_MZ: " << APFEL::AlphaQCD(APFEL::GetZMass()) << endl;
        infodata << scientific;
        infodata << "AlphaS_OrderQCD: " << APFEL::GetPerturbativeOrder() << endl;
        infodata << "AlphaS_Type: ipol" << endl;
        infodata << "AlphaS_Qs: [";

        for (int i = 0; i < nQ2; i++)
            infodata << sqrt(Q2_array[i]) << ((i == nQ2-1) ? "]\n" : ", ");

        infodata << "AlphaS_Vals: [";
        for (int i = 0; i < nQ2; i++)
            infodata << APFEL::AlphaQCD(sqrt(Q2_array[i])) << ((i == nQ2-1) ? "]\n" : ", ");

        infodata.close();
    }
  
    return 0;
}