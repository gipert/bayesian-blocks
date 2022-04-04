// Copyright (c) 2019 Luigi Pertoldi
//
// Permission is hereby granted, free of charge, to any person obtaining a copy of
// this software and associated documentation files (the "Software"), to deal in
// the Software without restriction, including without limitation the rights to
// use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
// the Software, and to permit persons to whom the Software is furnished to do so,
// subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
// FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
// COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
// IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
// CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#include "TH1.h"
#include "bayesian_blocks.hpp"

#ifndef _BAYESIAN_BLOCKS_ROOT_HH
#define _BAYESIAN_BLOCKS_ROOT_HH

namespace BayesianBlocks {

    // rebin a ROOT histogram
    TH1* rebin(TH1* h_in, const double p = 0.01,
               bool counter = false, bool benchmark = false);
}

namespace BayesianBlocks {

    TH1* rebin(TH1* h_in, const double p, bool counter, bool benchmark) {

        // define variables
        const auto Nb = h_in->GetNbinsX();
        bb::data_array x;
        bb::weights_array weights;
        bb::array edges;

        // find first non-empty bin
        int i_first = 1;
        for (int i = 1; i < Nb; ++i ) {
            if (h_in->GetBinContent(i) != 0) {
                edges  .push_back(h_in->GetBinCenter(i));
                x      .push_back(h_in->GetBinCenter(i));
                weights.push_back(h_in->GetBinContent(i));
                i_first = i;
                break;
            }
        }

        // fill arrays, skip empty bins
        for (int i = i_first+1; i < Nb; ++i ) {
            auto c = h_in->GetBinContent(i);
            if (floor(c) != ceil(c)) {
                throw std::domain_error("ERROR: non-integer bin contents detected in input histogram");
            }
            if (c == 0) continue;
            x      .push_back(h_in->GetBinCenter(i));
            weights.push_back(c);
        }

        auto result = BayesianBlocks::blocks(x, weights, p, counter, benchmark);

        auto h_out = new TH1D(
            (std::string(h_in->GetName()) + "_b").c_str(), h_in->GetTitle(),
            result.size()-1, &result[0]
        );

        for (int b = 1; b < h_in->GetNbinsX(); ++b) {
            auto c = h_in->GetBinContent(b);
            auto bin = h_out->FindBin(h_in->GetBinCenter(b));
            h_out->SetBinContent(bin, h_out->GetBinContent(bin) + c);
        }
        h_out->SetBinContent(0, h_in->GetBinContent(0)); // underflow bin
        h_out->SetBinContent(result.size(), h_in->GetBinContent(h_in->GetNbinsX())); // overflow bin
        h_out->Scale(1, "width"); // make it a density by dividing contents and errors by bin widths

        return h_out;
    }
}

#endif
