#include "NuXSSplMkr/spline_maker.h"

namespace nuxssplmkr {

SplineMaker::SplineMaker(){ 

}

void SplineMaker::MakeSpline(string infile, string outfile) {
    std::ifstream gridfile(infile);
    std::string line;

    // The first line contains the log Q2 values
    std::getline(gridfile, line);
    std::istringstream logQ2_stream(line);
    std::vector<double> logQ2_values;
    double _logQ2;
    while (logQ2_stream >> _logQ2) {logQ2_values.push_back(_logQ2);}

    // The second line contains the log x values
    std::getline(gridfile, line);
    std::istringstream logx_stream(line);
    std::vector<double> logx_values;
    double _logx;
    while (logx_stream >> _logx) {logx_values.push_back(_logx);}

    unsigned int nQ2 = logQ2_values.size();
    unsigned int nx = logx_values.size();
    std::deque<std::pair<double,std::array<unsigned int, 2>>> data;
    photospline::ndsparse spline_data(nx*nQ2, 2);
    std::vector<double> weights(nx*nQ2, 1.);

    // Construct the knots
    std::vector<double> Q2_knots;
    std::vector<double> x_knots;

    // Get bounds from the grids
    double log_Q2_min = logQ2_values.at(0);
    double log_Q2_max = logQ2_values.at(nQ2-1);
    double log_x_min = logx_values.at(0);
    double log_x_max = logx_values.at(nx-1);

    const uint32_t dim = 2;
    unsigned int Nknots_Q2 = nx + 2;
    unsigned int Nknots_x = nQ2 + 2;
    const double d_log_Q2_knot= std::abs(log_Q2_max - log_Q2_min) / (Nknots_Q2 - 1);
    const double d_log_x_knot = std::abs(log_x_max - log_x_min ) / (Nknots_x - 1);

    // Need to make some assumptions here
    // ### Testing new knots ###
    for ( double log_Q2 = log_Q2_min - 1.0; log_Q2 <= log_Q2_min; log_Q2 += 0.1 ) {
        double knot = log_Q2;
        Q2_knots.push_back(knot);
    }

    for ( double log_Q2 = log_Q2_min + d_log_Q2_knot; log_Q2 <= log_Q2_max + d_log_Q2_knot; log_Q2 += d_log_Q2_knot ) {
        double knot = log_Q2;
        Q2_knots.push_back(knot);
    }

    for ( double log_x = log_x_min - d_log_x_knot; log_x <= std::log10(0.5)-d_log_x_knot; log_x += d_log_x_knot ) {
        double knot = log_x;
        x_knots.push_back(knot);
    }

    // for ( double log_x = -1.0; log_x <= 0.25; log_x += 0.01 ) {
    //     double knot = log_x;
    //     x_knots.push_back(knot);
    // }

    for ( double x = 0.5; x < 0.9; x += 0.025 ) {
        double knot = std::log10(x);
        x_knots.push_back(knot);
    }

    for ( double x = 0.9; x < 1.2; x += 0.01 ) {
        double knot = std::log10(x);
        x_knots.push_back(knot);
    }

    for ( double x = 1.2; x <= 2.0; x += 0.05 ) {
        double knot = std::log10(x);
        x_knots.push_back(knot);
    }
    // ####

    double smooth = 1e-15;
    std::vector<uint32_t> orders(dim, 2);

    for (unsigned int i = 0; i < nQ2; i++){
        std::getline(gridfile, line);
        std::istringstream linestream(line);
        std::string val;
        unsigned int j = 0;
        while (std::getline(linestream, val, ',')) {
            data.push_back(std::make_pair(std::stod(val), std::array<unsigned int, 2>{i, j}));
            // std::cout << logQ2_values.at(i) << "," << logx_values.at(j) << ": " << val << std::endl;
            j += 1;
        }
    }
    for(auto& entry : data) {
        spline_data.insertEntry(entry.first, &entry.second[0]);
    }

    photospline::splinetable<> spline;
    spline.fit(spline_data, weights, std::vector<std::vector<double>>{logQ2_values, logx_values}, orders, 
               {Q2_knots, x_knots}, {smooth, smooth}, {2, 2});
    if(std::isnan(*spline.get_coefficients())){
        std::cerr << "Spline fit has failed!" << std::endl;
    }
    spline.write_fits(outfile);
}

}