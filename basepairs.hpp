#ifndef INCL_BASEPAIRS
#define INCL_BASEPAIRS

#include <string>
#include <fstream>

#define ALPH_SIZE 4

typedef unsigned char base_t;

std::string reverse_complement(const std::string &seq);
    
base_t char_to_base(char b);

char base_to_char(base_t i);

void parse_fasta(std::ifstream &fasta_in, 
                 std::string &fwd_bases, 
                 std::string &rev_bases);

#endif
