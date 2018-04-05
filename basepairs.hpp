#ifndef INCL_BASEPAIRS
#define INCL_BASEPAIRS

#include <string>
#include <fstream>
#include <vector>

#define ALPH_SIZE 4

typedef unsigned char Base;

std::string reverse_complement(const std::string &seq);

std::vector<Base> seq_to_bases(const std::string &seq);
    
Base char_to_base(char b);

char Baseo_char(Base i);

void parse_fasta(std::ifstream &fasta_in, 
                 std::string &fwd_bases, 
                 std::string &rev_bases,
                 bool terminate=false);

#endif
