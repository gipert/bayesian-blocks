#include <iostream>
#include <fstream>
#include <vector>

#include "../bayesian_blocks.hh"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"

int main() {

    std::ifstream fin("test.dat");
    double a;
    std::vector<double> v;
    while (fin >> a) v.emplace_back(a);

    auto r = BayesianBlocks::blocks(v, 0.01, false, true);
    std::cout << "Bin edges: ";
    for (auto& i : r) std::cout << i << " "; std::cout << std::endl;

    std::vector<double> w = {1, 2, 3, 3, 2, 7, 6, 5, 9, 6, 5, 5, 6, 4, 7, 3, 5, 4, 3, 6, 5, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};
    auto rr = BayesianBlocks::blocks(w, 0.01, false, true);
    std::cout << "Bin edges: ";
    for (auto& i : rr) std::cout << i << " "; std::cout << std::endl;

    // testing with ROOT histogram
    TFile f("test.root");
    auto h = dynamic_cast<TH1D*>(f.Get("hmodel_ch1"));
    h = dynamic_cast<TH1D*>(h->Rebin(2)); // or it takes too long...
    std::cout << "Rebinning..." << std::endl;
    auto hr = dynamic_cast<TH1D*>(BayesianBlocks::rebin(h, 0.1, true, true));
    // BayesianBlocks::rebin(h, 0.01, true);
    std::cout << "\nDone!\n";

    TCanvas c;
    h->GetYaxis()->SetRangeUser(1E-5, 1E4);
    h->Draw("HIST");
    hr->Draw("HIST SAME");
    c.SetLogy();
    c.SaveAs("test.pdf");
}
