#include <iostream>
#include <string>
#include <unordered_map>
#include <math.h>
#include <unistd.h>

//#include "kmer_model.cpp"
#include "sdsl_fmi.hpp"
#include "base_fmi.hpp"
#include "timer.hpp"

int main(int argc, char** argv) {

    std::string ref_fname(argv[1]);
    std::string out_fname(argv[2]);
    //std::string suffix = "idx";


    std::string combined_bases, rev_bases;

    std::ifstream ref_file(ref_fname);
    parse_fasta(ref_file, combined_bases, rev_bases, false);

    combined_bases.insert(combined_bases.end(), 
                          rev_bases.begin(), 
                          rev_bases.end());
    combined_bases = combined_bases + "$";

    std::cerr << "Building forward FMI\n";
    SdslFMI fmi;
    fmi.construct(combined_bases);
    //BaseFMI fmi(combined_bases, 128);

    std::cerr << "Saving combined FMI\n";
    fmi.save(out_fname);

    //std::cerr << "Building reverse FMI\n";
    //SdslFMI rev_fmi;
    //rev_fmi.construct(rev_bases);
    //BaseFMI rev_fmi(rev_bases, 128);

    //std::cerr << "Saving reverse FMI\n";
    //rev_fmi.save(out_prefix + "revFM." + suffix);

    return 0;
}
