#include <iostream>
#include <string>
#include <unordered_map>
#include <math.h>
#include <unistd.h>

#include "kmer_model.cpp"
#include "base_fmi.hpp"
#include "timer.h"

int main(int argc, char** argv) {

    std::string ref_fname(argv[1]);
    std::string out_prefix(argv[2]);

    std::ifstream ref_file(ref_fname);

    std::string fwd_bases, rev_bases;
    parse_fasta(ref_file, fwd_bases, rev_bases);

    std::cerr << "Building forward FMI\n";
    BaseFMI fwd_fmi(fwd_bases, (unsigned int) 200);

    std::cerr << "Saving forward FMI\n";
    fwd_fmi.save(out_prefix + "fwdFM.txt");

    std::cerr << "Building reverse FMI\n";
    BaseFMI rev_fmi(rev_bases, (unsigned int) 200);

    std::cerr << "Saving reverse FMI\n";
    rev_fmi.save(out_prefix + "revFM.txt");

    return 0;
}
