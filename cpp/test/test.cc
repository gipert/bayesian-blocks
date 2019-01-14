#include <iostream>
#include <vector>

#include "bayesian_blocks.hh"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"

int main() {

    // testing with two simple vectors
    std::vector<double> v = {1, 3, 2, 5, 4, 3, 8, 6, 4, 5, 7, 4, 3, 6, 5, 2, 6, 5, 4};
    std::vector<double> w = {1, 3, 2, 5, 4, 3, 8, 6, 4, 5, 7, 4, 3, 6, 5, 2, 6, 5, 4};

    auto r = BayesianBlocks::blocks(v, w, 0.01);
    std::cout << "Bin edges: ";
    for (auto& i : r) std::cout << i << " "; std::cout << std::endl;

    // testing with ROOT histogram
    TFile f("test.root");
    auto h = dynamic_cast<TH1D*>(f.Get("hmodel_ch1"));
    h = dynamic_cast<TH1D*>(h->Rebin(2)); // or it takes too long...
    std::cout << "Rebinning..." << std::endl;
    auto hr = dynamic_cast<TH1D*>(BayesianBlocks::rebin(h, 0.01, true));
    std::cout << "\nDone!\n";

    TCanvas c;
    h->GetYaxis()->SetRangeUser(1E-5, 1E4);
    h->Draw("HIST");
    hr->Draw("SAME HIST");
    c.SetLogy();
    c.SaveAs("test.pdf");
}
