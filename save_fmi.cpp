#include <iostream>
#include <string>
#include <unordered_map>
#include <math.h>
#include <unistd.h>

#include "kmer_model.cpp"
#include "base_fmi.hpp"
#include "timer.h"

std::string reverse_complement(const std::string &seq) {
    std::string rev(seq.size(), 'N');
    for (unsigned int i = 0; i < seq.size(); i++) {
        char c = 'N';
        switch(seq[i]) {
            case 'A':
            case 'a':
            c = 'T';
            break;
            case 'T':
            case 't':
            c = 'A';
            break;
            case 'G':
            case 'g':
            c = 'C';
            break;
            case 'C':
            case 'c':
            c = 'G';
            break;
        }
        rev[rev.size()-i-1] = c;
    }

    return rev;
}

void parse_fasta(std::ifstream &fasta_in, 
                 std::string &fwd_bases, 
                 std::string &rev_bases) {

    //For parsing the file
    std::string line;
    
    getline(fasta_in, line); //read past header

    while (getline(fasta_in, line)) {
        fwd_bases += line;
    }

    rev_bases = reverse_complement(fwd_bases);

    fwd_bases += "$";
    rev_bases += "$";
}

int main(int argc, char** argv) {

    std::string ref_fname(argv[1]);
    std::string out_prefix(argv[2]);

    std::ifstream ref_file(ref_fname);

    std::string fwd_bases, rev_bases;
    parse_fasta(ref_file, fwd_bases, rev_bases);

    std::cerr << "Building forward FMI\n";
    BaseFMI fwd_fmi(fwd_bases, 200);

    std::cerr << "Saving forward FMI\n";
    fwd_fmi.save(out_prefix + "fwdFM.txt");

    std::cerr << "Building reverse FMI\n";
    BaseFMI rev_fmi(rev_bases, 200);

    std::cerr << "Saving reverse FMI\n";
    rev_fmi.save(out_prefix + "revFM.txt");

    return 0;
}
