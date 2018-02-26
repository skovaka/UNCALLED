#include <iostream>
#include <string>
#include <list>
#include "base_fmi.hpp"

int main(int argc, char** argv) {
    std::string ref = "GGCCTACTACATGCATGCAGTCCAAGTGGGCCCAATTCACTAACACTGATTACAACACTACGA$";
    std::string read = "GATTACA";

    int t = std::atoi(argv[1]);

    BaseFMI fmi(ref, t);

    Range r = fmi.get_full_range(read.back());

    for (int s = r.start_; s <= r.end_; s++) {
        std::cout << (fmi.suffix_ar_[s]) << "(" << ref[fmi.suffix_ar_[s]] << ") ";
    }
    std::cout << "\n";

    for (size_t i = read.size()-2; i > 0; i--) {

        r = fmi.get_neighbor(r, read[i]);

        for (int s = r.start_; s <= r.end_; s++) {
            std::cout << (fmi.suffix_ar_[s]) << "(" << s << "," << ref[fmi.suffix_ar_[s]] << ") ";
        }
        std::cout << "\n";
    }


}
