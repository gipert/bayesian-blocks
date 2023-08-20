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

// parse items in the format file.root:object
std::pair<std::string, std::string> get_file_obj(std::string expr) {
    std::string filename;
    std::string objname = "";
    if (expr.find(':') != std::string::npos) {
        filename = expr.substr(0, expr.find_first_of(':'));
        objname = expr.substr(expr.find_first_of(':')+1, std::string::npos);
    }
    else filename = expr;

    return std::pair<std::string, std::string>(filename, objname);
}

// splits filename from path
std::pair<std::string, std::string> split_dir_file(std::string expr) {
    std::string directory, file;
    if (expr.back() == '/') expr.pop_back();
    if (expr.find('/') != std::string::npos) {
        directory = expr.substr(0, expr.find_last_of('/'));
        file = expr.substr(expr.find_last_of('/')+1, std::string::npos);
    }
    else file = expr;

    return std::pair<std::string, std::string>(directory, file);
}

int main(int argc, char** argv) {

    std::string progname(argv[0]);

    auto usage = [&]() {
        std::cerr << "USAGE: " << progname << " [-v|--verbose] [-h|--help] [--p0 <val> (default 0.01)] file|file:obj [file2...]\n";
    };

    const char* const short_opts = "o:usvh";
    const option long_opts[] = {
        { "p0",      required_argument, nullptr, 1   },
        { "output",  required_argument, nullptr, 'o' },
        { "uniform", no_argument,       nullptr, 'u' },
        { "samedir", no_argument,       nullptr, 's' },
        { "verbose", no_argument,       nullptr, 'v' },
        { "help",    no_argument,       nullptr, 'h' },
        { nullptr,   no_argument,       nullptr, 0   }
    };

    // defaults
    double p0 = 0.01;
    std::string outfile = "";
    bool uniform = false;
    bool samedir = false;

    int opt = 0;
    while ((opt = getopt_long(argc, argv, short_opts, long_opts, nullptr)) != -1) {
        switch (opt) {
            case 1:
                p0 = std::stod(optarg);
                break;
            case 'o':
                outfile = std::string(optarg);
                break;
            case 'u':
                uniform = true;
                break;
            case 's':
                samedir = true;
                break;
            case 'v':
                level = debug;
                break;
            case '?': // Unrecognized option
                glog(error) << "ERROR: unrecognized option!" << std::endl;
                return 1;
            case 'h': // -h or --help
            default:
                usage();
                return 1;
        }
    }

    if (p0 != 0.01) glog(debug) << "p0 set to " << p0 << std::endl;
    if (!outfile.empty()) glog(debug) << "custom output file set to " << outfile << std::endl;

    // extra arguments
    std::vector<std::string> args;
    for(; optind < argc; optind++){
        args.emplace_back(argv[optind]);
    }
    if (args.empty()) {usage(); return 1;}

    TArrayD saved_edges;

    // iterate over files
    // if no output filename is specified make an output file for each file
    for (auto& f : args) {
        glog(debug) << f << std::endl;

        std::vector<TObject*> hists;
        auto file_obj = get_file_obj(f);

        // open file
        TFile _tmp(file_obj.first.c_str(), "read");

        // if no object is specified, try processing all histograms in the file
        // will be saved in the same output file
        if (file_obj.second.empty()) {
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
        }
        // otherwise, process that object only
        else {
            auto h = dynamic_cast<TH1*>(_tmp.Get(file_obj.second.c_str()));
            if (!h) throw std::runtime_error("Could not read object '" + file_obj.second + "' in file as histogram");

            TString tmp_file_obj(file_obj.second);
            tmp_file_obj.ReplaceAll('/', '_');
            h->SetName(tmp_file_obj);

            glog(debug) << " ├─ " << h->GetName();

            TH1* hr;
            // use edges from first histogram
            if (saved_edges.GetSize() > 0 and uniform == true) {
                hr = h->Rebin(saved_edges.GetSize()-1, tmp_file_obj, saved_edges.GetArray());
                hr->Scale(1, "width");
                std::cout << " (using cached bin edges)";
            }
            // compute edges
            else {
                bool found = false;
                for (int b = 1; b < h->GetNbinsX(); ++b) {
                    auto c = h->GetBinContent(b);
                    if (floor(c) != ceil(c)) {
                        if (!found) std::cerr << " \033[93m✘\033[0m non-integer bin contents detected, they will be rounded.";
                        h->SetBinContent(b, std::round(c));
                        found = true;
                    }
                }

                try {
                    hr = dynamic_cast<TH1D*>(BayesianBlocks::rebin(h, p0, false, false));
                    saved_edges = *hr->GetXaxis()->GetXbins();
                }
                catch(const std::exception& e) {
                    if (level <= debug) std::cerr << " \033[91m✘\033[0m " << e.what() << ". Skipping.\n";
                    continue;
                }
            }
            hists.push_back(hr);
            if (level <= debug) std::cout << " \033[92m✔\033[0m\n";
        }
        glog(debug) << " └─ done\n";

        if (!hists.empty()) {
            std::string outname;
            if (!outfile.empty())
                outname = outfile;
            else {
                outname = "bb-" + (split_dir_file(file_obj.first)).second;
                if (samedir) outname = (split_dir_file(file_obj.first)).first + "/" + outname;
            }
            TFile fout(outname.c_str(), "update");
            for (auto& i : hists) {
                auto name = std::string(i->GetName());
                if (*(name.end()-1) == 'b' and *(name.end()-2) == '_') {
                    name.erase(name.end()-2, name.end());
                }
                i->Write(name.c_str());
            }
            glog(info) << outname << " written\n";
        }
    }
    return 0;
}
