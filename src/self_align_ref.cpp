#include <iostream>
#include <string>
#include <list>
#include <cstdlib>
#include "range.hpp"
#include "bwa_fmi.hpp"
#include "timer.hpp"
#include "self_align_ref.hpp"

std::vector< std::vector<u64> > self_align(const std::string &bwa_prefix,
                                         const std::string fasta_fname,
                                         u32 sample_dist) {

    BwaFMI fmi(bwa_prefix);

    std::vector< std::vector<u8> > seqs;
    std::ifstream fasta_in(fasta_fname);
    std::string fasta_line;
    while (getline(fasta_in, fasta_line)) {
        if (fasta_line[0] == '>') {
            seqs.push_back(std::vector<u8>());
        } else {
            for (char c : fasta_line) {
                seqs.back().push_back(BASE_COMP_B[BASE_BYTES[(u8)c]]);
            }
        }
    }

    srand(0);
    
    std::vector< std::vector<u64> > ret;

    for (auto bases : seqs) {
        for (u64 i = 0; i < bases.size(); i++) {
            if (rand() % sample_dist != 0) {
                continue;
            }

            ret.push_back(std::vector<u64>());

            Range r = fmi.get_full_range(bases[i]);
            u64 j = i+1;
            for (; j < bases.size() && r.length() > 1; j++) {
                ret.back().push_back(r.length());
                r = fmi.get_neighbor(r, bases[j]);
            }
            ret.back().push_back(r.length());
        }
    }

    return ret;
}

