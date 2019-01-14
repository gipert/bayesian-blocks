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

#include "bayesian_blocks.hh"

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <stdexcept>

#include "TH1.h"

namespace BayesianBlocks {

    bb::vec blocks(bb::vec data, bb::vec weights, const double p,
                   bool counter, bool benchmark) {

        auto start = bb::clock::now();

        // sanity checks
        if (data.size() != weights.size()) {
            throw std::domain_error("data and weights vectors are of different sizes");
        }

        if (std::find_if(weights.begin(), weights.end(),
                         [](double& v) { return v <= 0; }) != weights.end()) {
            throw std::domain_error("Invalid weights found in input");
        }

        const auto N = data.size();

        // sort and copy data
        std::vector<std::pair<double&, double&>> hist; hist.reserve(N);
        for (std::size_t i = 0; i < N; ++i) hist.emplace_back(data[i], weights[i]);
        std::sort(hist.begin(), hist.end(), [](bb::pair a, bb::pair b) { return a.first < b.first; });

        // build up array with all possible bin edges
        bb::vec edges(N+1);
        for (std::size_t i = 0; i <= N; ++i) {
            edges[i] =
                ( i == 0 ? edges[0] = data[0] :
                ( i == N ? edges[N] = data[N-1] :
                (data[i]+data[i+1])/2.));
        }

        // let's use here Cash statistics and calibrated prior on number of change points
        auto cash = [](int N_k, double T_k) { return N_k * std::log(N_k/T_k); };
        auto ncp_prior = std::log(73.53 * p * std::pow(N, -0.478)) - 4;

        // arrays to store results
        bb::vec last(N);
        bb::vec best(N);

        auto init_time = bb::duration_cast<bb::us>(bb::clock::now() - start).count();
        start = bb::clock::now();

        // do the actual recursive computation
        for (std::size_t k = 0; k < N; ++k ) {
            bb::vec A(k+1);
            for (std::size_t r = 0; r == 0 or r <= k; ++r) {
                A[r] = cash(
                            std::accumulate(weights.begin()+r, weights.begin()+k+1, 0),
                            edges[k+1] - edges[r]
                       ) + ncp_prior + (r == 0 ? 0 : best[r-1]);
            }
            last[k] = std::distance(A.begin(), std::max_element(A.begin(), A.end()));
            best[k] = *(std::max_element(A.begin(), A.end()));

            if (counter) std::cout << '\r' << k << '/' << N << std::flush;
        }
        if (counter) std::cout << std::endl;

        auto loop_time = bb::duration_cast<bb::us>(bb::clock::now() - start).count();
        start = bb::clock::now();

        // iteratively find the change points
        std::vector<int> cp;
        for (auto i = N; i != 0; i = last[i-1]) cp.push_back(i); cp.push_back(0);

        std::reverse(cp.begin(), cp.end());
        bb::vec result(cp.size(), 0);
        std::transform(cp.begin(), cp.end(), result.begin(),
                       [edges](size_t pos) { return edges[pos]; });

        auto end_time = bb::duration_cast<bb::us>(bb::clock::now() - start).count();

        if (benchmark) {
            std::cout << "init: ";
            init_time > 1000 ?
                std::cout << init_time/1.E3 << " s" :
                std::cout << init_time      << " us";
            std::cout << std::endl;
            std::cout << "loop: ";
            loop_time > 1000 ?
                std::cout << loop_time/1.E3 << " s" :
                std::cout << loop_time      << " us";
            std::cout << std::endl;
            std::cout << "end: ";
            end_time > 1000 ?
                std::cout << end_time/1.E3 << " s" :
                std::cout << end_time      << " us";
            std::cout << std::endl;
        }

        return result;
    }

    TH1* rebin(TH1* h_in, const double p, bool counter) {

        // define variables
        const auto Nb = h_in->GetNbinsX();
        bb::vec x;
        bb::vec weights;
        bb::vec edges;

        // find first non-empty bin
        int i_first = 1;
        for (int i = 1; i <= Nb; ++i ) {
            if (h_in->GetBinContent(i) != 0) {
                edges  .push_back(h_in->GetBinCenter(i));
                x      .push_back(h_in->GetBinCenter(i));
                weights.push_back(h_in->GetBinContent(i));
                i_first = i;
                break;
            }
        }

        // fill arrays, skip empty bins
        for (int i = i_first+1; i <= Nb; ++i ) {
            if (h_in->GetBinContent(i) == 0) continue;
            x      .push_back(h_in->GetBinCenter(i));
            weights.push_back(h_in->GetBinContent(i));
        }

        auto result = BayesianBlocks::blocks(x, weights, p, counter);

        auto h_out = dynamic_cast<TH1*>(
            h_in->Rebin(result.size()-1,
                        (std::string(h_in->GetName()) + "_b").c_str(),
                        &result[0]
            )
        );
        h_out->Scale(1, "width");

        return h_out;
    }
}

// vim: foldmethod=syntax
