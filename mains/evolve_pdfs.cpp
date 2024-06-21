#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cmath>
#include <boost/filesystem.hpp>
#include "APFEL/APFEL.h"

using namespace std;

int main(int argc, char* argv[]){
    if (argc < 5) {
        std::cout << "Not enough inputs." << std::endl;
        return 0;
    }
    std::string inpdf = argv[1];  // input pdfset name
    std::string outpdf = argv[2]; // output pdfset name
    std::string replica = argv[3];
    int rep = std::stoi(replica); // replica number
    int pto = std::stoi(argv[4]); // perturbative order
    bool enable_small_x = !!(std::stoi(argv[5])); // use small-x resummation?

    auto replica_string = std::string(4 - std::min(4, static_cast<int>(replica.length()) ), '0') + replica;
    boost::filesystem::path out_folder = outpdf+"/";
    if (!boost::filesystem::exists(out_folder)) {
        boost::filesystem::create_directories(out_folder);
    }
    std::ofstream outfile;
    outfile.open(outpdf + "/" + outpdf + "_" + replica_string + ".dat");

    if (replica == "0") {
        outfile << "PdfType: central" << std::endl;
    } else {
        outfile << "PdfType: replica" << std::endl;
    }
    outfile << "Format: lhagrid1" << std::endl << "---" << std::endl;

    double Q0 = 1.295;
    if (enable_small_x) {
        Q0 = std::sqrt(3.5); // TODO: be smarter about this? -PW
    }

    APFEL::SetPDFSet(inpdf);
    APFEL::SetReplica(0);
    APFEL::SetTheory("QCD");
    APFEL::SetQLimits(Q0, 1e5);
    APFEL::SetPerturbativeOrder(pto);
    APFEL::SetPDFEvolution("exactalpha");
    APFEL::SetFastEvolution(false);
    APFEL::SetAlphaEvolution("exact");
    APFEL::SetVFNS();
    APFEL::SetMaxFlavourPDFs(6);
    APFEL::SetMaxFlavourAlpha(6);
    APFEL::SetNumberOfGrids(3);
    APFEL::SetZMass(91.1870);
    APFEL::SetGridParameters(1, 130, 3, 1e-9);
    APFEL::SetGridParameters(2, 60, 5, 1e-1);
    APFEL::SetGridParameters(3, 20, 5, 8e-1);
    APFEL::SetPoleMasses(1.3, 4.75, 172.0);
    APFEL::SetAlphaQCDRef(0.118, 91.1870);

    int nQ2 = 100;
    int nlogx = 100;
    int nlinx = 50;
    int nx = nlogx + nlinx;

    double logQ2_min = std::log10(Q0*Q0);
    double logQ2_max = 10;
    double logx_min = -9;
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

    double pdf[nQ2][nx][13];
    for (int i = 0; i < nQ2; i++) {
        APFEL::EvolveAPFEL(Q0, std::sqrt(Q2_array[i]));
        // std::cout << Q2_array[i] << std::endl;

        for (int j = 0; j < nx; j++) {
            double x = x_array[j];
            for (int p = -6; p < 7; p++) {
                double xf = APFEL::xPDF(p, x);
                pdf[i][j][p+6] = xf;
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
                outfile << pdf[i][j][p+6] << " ";
            }
            outfile << std::endl;
        }
    }
    outfile << "---";
    outfile.close();

    if (!boost::filesystem::exists(outpdf + "/" + outpdf + ".info")) {
        // The following code is from the nnpdf collab code!
        // LHAPDF6 HEADER
        std::ofstream infodata;
        infodata.open(outpdf + "/" + outpdf + ".info");
        int nrep = 1; // TODO!
        infodata << "SetDesc: \"test\"" << endl;
        infodata << "SetIndex: " << endl;
        infodata << "Authors: Philip Weigel" << endl;
        infodata << "Reference: N/A" << endl;
        infodata << "Format: lhagrid1" << endl;
        infodata << "DataVersion: 1" << endl;
        infodata << "NumMembers: " << nrep << endl;
        infodata << "Particle: 2212" << endl;
        infodata << "Flavors: [";
      for (int p = -6; p < 7; p++) {
            if (p == 0) {
                infodata << "21, ";
            } else if (p == 6) {
                infodata << "6";
            } else{
                infodata << p << ", ";
            }
        }
        infodata << "]" << endl;
        infodata << "OrderQCD: " << APFEL::GetPerturbativeOrder() << endl;
        infodata << "FlavorScheme: variable" << endl;
        infodata << "NumFlavors: " << std::max(APFEL::GetMaxFlavourPDFs(), APFEL::GetMaxFlavourAlpha()) << endl;
        infodata << "ErrorType: replicas" << endl;

        infodata.precision(7);
        infodata << scientific;
        infodata << "XMin: "<< x_array[0] << endl;
        infodata << "XMax: "<< x_array[nx-1] << endl;
        infodata << "QMin: "<< Q0 << endl;
        infodata << "QMax: "<< sqrt(pow(10, logQ2_max)) << endl;
        infodata << "MZ: "  << APFEL::GetZMass() << endl;
        infodata << "MUp: 0\nMDown: 0\nMStrange: 0" << std::endl;
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

        infodata << "AlphaS_Lambda4: 0.342207" << endl;
        infodata << "AlphaS_Lambda5: 0.239" << endl;
        infodata.close();
    }
  
    return 0;
}