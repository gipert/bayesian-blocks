#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>

#include "../bayesian_blocks.hh"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"

bool equal(std::vector<double> a, std::vector<double> b, double eps);

int main() {

    int exit_status = 0;
    std::ifstream fin("test.dat");
    double a;
    std::vector<double> v;
    while (fin >> a) v.emplace_back(a);

    std::vector<double> exp = {
        -3.48528, -1.87114, -1.36282, -0.677218,
        0.659105, 1.39771, 4.06582, 5.60912,
        6.17286, 7.76634, 9.91696
    };
    auto r = BayesianBlocks::blocks(v, 0.01, false, true);
    if (!equal(r, exp, 1E-04)) {
        std::cerr << "ERROR: test 1 failed.\n";
        std::cerr << "Expected: ";
        for (auto& val : exp) std::cerr << val << ' '; std::cerr << std::endl;
        std::cerr << "Got: ";
        for (auto& val : r) std::cerr << val << ' '; std::cerr << std::endl;
        exit_status = 1;
    }

    // testing with ROOT histogram
    TH1D h("hist", "hist", 100, -3, 10);
    for (auto& val : v) h.Fill(val);
    auto hr = dynamic_cast<TH1D*>(BayesianBlocks::rebin(&h, 0.01, false, true));
    auto edges = hr->GetXaxis()->GetXbins()->GetArray();
    std::vector<double> w(edges, edges + hr->GetNbinsX()+1);

    exp = {
        -2.935, -1.7, -1.18, -0.66, 0.64,
        1.42, 4.02, 5.58, 6.23, 7.79, 9.805
    };
    if (!equal(w, exp, 1E-03)) {
        std::cerr << "ERROR: test 2 failed:" << std::endl;
        std::cerr << "Expected: ";
        for (auto& val : exp) std::cerr << val << ' '; std::cerr << std::endl;
        std::cerr << "Got: ";
        for (auto& val : w) std::cerr << val << ' '; std::cerr << std::endl;
        exit_status = 2;
    }

    return exit_status;
}

bool equal(std::vector<double> a, std::vector<double> b, double eps) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (std::abs(a[i] - b[i]) > eps) return false;
    }
    return true;
}
