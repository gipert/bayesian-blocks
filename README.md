# Bayesian Blocks

Header-only C++ implementation of the Bayesian Blocks argorithm. Plug and play!

```cpp
#include "baeysian_blocks.hpp"
#include "baeysian_blocks_root.hpp"  // if you want to rebin ROOT CERN histograms

int main() {

    // "data" is a a std::vector containing the observed values
    // [optional] "weights" can be used to specify the weight (how many times an observation occurs)
    auto change_points = BayesianBlocks::blocks(data, /* optional */ weights);
    // "change_points" is a std::vector
 
    // "hist" is a ROOT TH1*
    auto h_rebin = dynamic_cast<TH1*>(BayesianBlocks::rebin(hist));
    // "h_rebin" is the rebinned histogram
    
    return 0;
}
```
Check out the `BayesianBlocks::blocks` and `BayesianBlocks::rebin` signatures for a more advanced usage. Have a look at [`test/run_test.cc`](https://github.com/gipert/bayesian-blocks/blob/master/test/run_test.cc) too.

Bayesian blocks algorithmScargle, J et al. (2012) [https://doi.org/10.1088/0004-637X/764/2/167] reference: *Scargle, J et al. (2012) [https://doi.org/10.1088/0004-637X/764/2/167]*

## The `bblocks` utility

The `bblocks` program rebins all the 1-dim histograms contained in a ROOT file.

Install:
```console
$ cd cpp
$ PREFIX=<install-prefix> make install
```

Run:
```console
$ bblocks --help
USAGE: bblocks [-v|--verbose] [-h|--help] [--p0 <val> (default 0.01)] file1 file2 ...
```

### Related
Julia language enthusiast? Checl out my Julia package: https://github.com/gipert/BayesianBlocks.jl
