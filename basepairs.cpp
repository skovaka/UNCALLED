
#include <string>
#include <fstream>
#include "basepairs.hpp"

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


    
base_t char_to_base(char b) {
    switch (b) {
        case 'A':
        case 'a':
        return 0;

        case 'C':
        case 'c':
        return 1;

        case 'G':
        case 'g':
        return 2;

        case 'T':
        case 't':
        return 3;
    }

    return -1;
}

char base_to_char(base_t i) {
    switch (i) {
        case 0:
        return 'A';
        
        case 1:
        return 'C';

        case 2:
        return 'G';

        case 3:
        return 'T';
    }

    return 'N';
}

std::vector<base_t> seq_to_bases(const std::string &seq) {
    std::vector<base_t> bases(seq.size());
    for (size_t i = 0; i < seq.size(); i++) {
        bases[i] = char_to_base(seq[i]);
    }
    return bases;
}

void parse_fasta(std::ifstream &fasta_in, 
                 std::string &fwd_bases, 
                 std::string &rev_bases,
                 bool terminate) {

    //For parsing the file
    std::string line;
    
    getline(fasta_in, line); //read past header

    while (getline(fasta_in, line)) {
        fwd_bases += line;
    }

    rev_bases = reverse_complement(fwd_bases);
    
    if (terminate) {
        fwd_bases += "$";
        rev_bases += "$";
    }
}
