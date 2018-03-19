#include <iostream>
#include <string>
#include <list>
#include "base_fmi.hpp"
#include "sdsl_fmi.hpp"
#include "fmi.hpp"
#include "timer.h"

void test(FMI &fmi, const std::vector<base_t> &bases) {
    Range r = fmi.get_full_range(bases.back());

    for (size_t i = bases.size()-2; i < bases.size(); i--) {
        //std::cout << (int) bases[i+1] << r.start_ << "-" << r.end_ << " " << r.length() << "\n";
        r = fmi.get_neighbor(r, bases[i]);
    }


    std::cout << "done\n";

    //for (size_t s = r.start_; s <= r.end_; s++) {
    //    std::cout << s << " " << fmi.sa(s) << "\n";
    //}
}

int main(int argc, char** argv) {
    std::string ref_fname(argv[1]), index_fname, fwd_str, rev_str;

    std::ifstream ref_file(ref_fname);
    parse_fasta(ref_file, fwd_str, rev_str, false);
    rev_str.erase();

    std::vector<base_t> bases = seq_to_bases(fwd_str);

    BaseFMI fmi;
    
    if (argc > 2) {
        index_fname = std::string(argv[2]);
        std::ifstream index_file(index_fname);
        fmi = BaseFMI(index_file, 128);
        //fmi = SdslFMI(index_fname);
    } else {
        fmi = BaseFMI(fwd_str + "$", 128);
        //fmi.construct(fwd_str);
    }
    
    Timer t;
    t.reset();

    std::list<Range> ranges;
    //Range r = fmi.get_full_range(bases.back());
    ranges.push_front(fmi.get_full_range(bases.back()));
    for (size_t i = bases.size()-2; i < bases.size() && ranges.front().length() > 1; i--) {
        std::cout << ranges.front().length() << "\n";
        ranges.push_front( fmi.get_neighbor(ranges.front(), bases[i]) );
    }
}
