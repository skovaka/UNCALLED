#include <iostream>
#include <string>
#include <list>
#include "base_fmi.hpp"
//#include "sdsl_fmi.hpp"
#include "fmi.hpp"
#include "timer.h"

void brute_align(FMI &fmi, const std::vector<base_t> &bases) {
    for (size_t j = bases.size()-1; j < bases.size(); j--) {
        Range r = fmi.get_full_range(bases[j]);
        size_t i = j-1;
        for (; i < bases.size() && r.length() > 1; i--) {
            //std::cout << r.length() << "\t";
            r = fmi.get_neighbor(r, bases[i]);
        }
        //std::cout << r.length() << "\n";
        std::cout << (j - i) << "\n";
    }
}

void dyn_align(FMI &fmi, const std::vector<base_t> &bases) {
    std::list<Range> ranges;
    ranges.push_front(fmi.get_full_range(bases.back()));

    for (size_t i = bases.size()-2; i < bases.size() && ranges.front().length() > 1; i--) {
        //std::cout << ranges.front().length() << "\t";
        ranges.push_front( fmi.get_neighbor(ranges.front(), bases[i]) );
    }
    //std::cout << ranges.front().length() << "\n";
    std::cout << ranges.size() << "\n";

    for (size_t i = bases.size()-2; i < bases.size(); i--) {
        ranges.pop_back();
        
        Range r = fmi.get_full_range(bases[i]), l;
        l = r.split_range(ranges.back());

        size_t j = i - 1;
        auto m = ranges.rbegin();
        for (; j < bases.size() && (l.is_valid() || r.is_valid()) && m != ranges.rend(); m++) {

            if (l.is_valid()) {
                m->start_ = l.start_;
                l = fmi.get_neighbor(l, bases[j]);
            }

            if (r.is_valid()) {
                m->end_ = r.end_;
                r = fmi.get_neighbor(r, bases[j]);
            }

            //std::cout << m->length() << "\t";

            j--;
        }

        //for (; m != ranges.rend(); m++) {
        //    std::cout << m->length() << "\t";
        //}
        //ONLY KEEP ONE
        //if (m != ranges.rend()) {
        //    std::cout << "x\t";
        //}

        j = i - ranges.size();

        while (ranges.front().length() > 1 && j < bases.size()) {
            r = fmi.get_neighbor(ranges.front(), bases[j]);
            //if (!r.is_valid()) {
            //    break;
            //}
            //std::cout << r.length() << "\t";
            ranges.push_front(r);
            j--;
        }
        //std::cout << "\n";
        std::cout << ranges.size() << "\n";
    }
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
    dyn_align(fmi, bases);
    std::cerr << t.lap() << "\n";
    //brute_align(fmi, bases);
    //std::cerr << t.lap() << "\n";


}
