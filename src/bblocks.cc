#include <getopt.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "../include/bayesian_blocks_root.hpp"

// ROOT
#include "TObject.h"
#include "TFile.h"
#include "TH1.h"
#include "TIterator.h"
#include "TKey.h"
#include "TClass.h"

enum log_level {debug, info, warning, error};

log_level level = info;

std::ofstream devnull("/dev/null");
std::ostream& glog(log_level lvl) {
    if (lvl == debug and level <= debug) {
        std::cout << "\033[36m[ Debug:\033[0m ";
        return std::cout;
    }
    if (lvl == info and level <= info) {
        std::cout << "\033[97;1m[ Info:\033[0m ";
        return std::cout;
    }
    if (lvl == warning and level <= warning) {
        std::cerr << "\033[33m[ Warning:\033[0m ";
        return std::cerr;
    }
    if (lvl == error and level <= error) {
        std::cerr << "\033[91m[ Error:\033[0m ";
        return std::cerr;
    }
    else {
        return devnull;
    }
};

int main(int argc, char** argv) {

    std::string progname(argv[0]);

    auto usage = [&]() {
        std::cerr << "USAGE: " << progname << " [-v|--verbose] [-h|--help] [--p0 <val> (default 0.01)] file1 file2 ...\n";
    };

    const char* const short_opts = ":vh";
    const option long_opts[] = {
        { "p0",      required_argument, nullptr, 1   },
        { "verbose", no_argument,       nullptr, 'v' },
        { "help",    no_argument,       nullptr, 'h' },
        { nullptr,   no_argument,       nullptr, 0   }
    };

    // default p0 value
    double p0 = 0.01;

    int opt = 0;
    while ((opt = getopt_long(argc, argv, short_opts, long_opts, nullptr)) != -1) {
        switch (opt) {
            case 1:
                p0 = std::stod(optarg);
                break;
            case 'v':
                level = debug;
                break;
            case 'h': // -h or --help
            case '?': // Unrecognized option
            default:
                usage();
                return 1;
        }
    }

    if (p0 != 0.01) glog(debug) << "p0 set to " << p0 << std::endl;

    // extra arguments
    std::vector<std::string> args;
    for(; optind < argc; optind++){
        args.emplace_back(argv[optind]);
    }
    if (args.empty()) {usage(); return 1;}

    // iterate over files
    for (auto& f : args) {
        glog(debug) << f << std::endl;

        std::vector<TObject*> hists;

        // open file
        TFile _tmp(f.c_str(), "read");
        TIter next(_tmp.GetListOfKeys());
        TKey* key;
        while ((key = dynamic_cast<TKey*>(next()))) {
            auto cl = TClass::GetClass(key->GetClassName());
            if (cl->InheritsFrom("TH1")) {
                auto h = dynamic_cast<TH1*>(key->ReadObj());
                if (h->GetDimension() == 1) {
                    glog(debug) << " ├─ " << h->GetName();

                    bool found = false;
                    for (int b = 1; b < h->GetNbinsX(); ++b) {
                        auto c = h->GetBinContent(b);
                        if (floor(c) != ceil(c)) {
                            if (!found) std::cerr << " \033[93m✘\033[0m non-integer bin contents detected, they will be rounded.";
                            h->SetBinContent(b, std::round(c));
                            found = true;
                        }
                    }


                    TH1* hr;
                    try {
                        hr = dynamic_cast<TH1D*>(BayesianBlocks::rebin(h, p0, false, false));
                    }
                    catch(const std::exception& e) {
                        if (level <= debug) std::cerr << " \033[91m✘\033[0m " << e.what() << ". Skipping.\n";
                        continue;
                    }
                    hists.push_back(hr);
                    if (level <= debug) std::cout << " \033[92m✔\033[0m\n";
                }
            }
        }
        glog(debug) << " └─ done\n";

        if (!hists.empty()) {
            TFile fout(("bb-" + f).c_str(), "recreate");
            for (auto& i : hists) {
                auto name = std::string(i->GetName());
                name.erase(name.end()-2, name.end());
                i->Write(name.c_str());
            }
        }

        glog(info) << "bb-" + f << " created\n";
    }
    return 0;
}
