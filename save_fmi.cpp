#include <iostream>
#include <string>
#include <unordered_map>
#include <math.h>
#include <unistd.h>

#include "kmer_model.cpp"
#include "sdsl_fmi.hpp"
#include "base_fmi.hpp"
#include "timer.h"

int main(int argc, char** argv) {

    std::string ref_fname(argv[1]);
    std::string out_prefix(argv[2]);
    std::string suffix = "idx";

    std::ifstream ref_file(ref_fname);

    std::string fwd_bases, rev_bases;
    parse_fasta(ref_file, fwd_bases, rev_bases, false);

    fwd_bases.erase();

    std::cerr << "Building forward FMI\n";
    //SdslFMI fwd_fmi;
    //fwd_fmi.construct(fwd_bases);
    ////BaseFMI fwd_fmi(fwd_bases, 128);

    //std::cerr << "Saving forward FMI\n";
    //fwd_fmi.save(out_prefix + "fwdFM." + suffix);

    std::cerr << "Building reverse FMI\n";
    SdslFMI rev_fmi;
    rev_fmi.construct(rev_bases);
    //BaseFMI rev_fmi(rev_bases, 128);

    std::cerr << "Saving reverse FMI\n";
    rev_fmi.save(out_prefix + "revFM." + suffix);

    return 0;
}
