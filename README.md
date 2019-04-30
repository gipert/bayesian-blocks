# Bayesian Blocks

[![Build Status](https://travis-ci.org/gipert/bayesian-blocks.svg?branch=master)](https://travis-ci.org/gipert/bayesian-blocks)

### The `bblocks` utility

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
