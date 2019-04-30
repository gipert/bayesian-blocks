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

#include <vector>
#include <chrono>

#include "TH1.h"

#ifndef _BAYESIAN_BLOCKS_HH
#define _BAYESIAN_BLOCKS_HH

namespace BayesianBlocks {

    // handy aliases
    namespace bb {
        // data containers
        using array         = std::vector<double>;
        using data_array    = std::vector<double>;
        using weights_array = std::vector<int>;
        using pair          = std::pair<double, int>;

        // time
        using clock = std::chrono::high_resolution_clock;
        using std::chrono::duration_cast;
        using us = std::chrono::microseconds;
    }

    // core utility
    bb::array blocks(bb::data_array data, bb::weights_array weights, const double p = 0.01,
                     bool counter = false, bool benchmark = false);

    bb::array blocks(bb::data_array data, const double p = 0.01,
                     bool counter = false, bool benchmark = false);

    // rebin a ROOT histogram
    TH1* rebin(TH1* h_in, const double p = 0.01,
               bool counter = false, bool benchmark = false);
}

#endif
