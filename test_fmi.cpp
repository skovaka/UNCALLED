#include <iostream>
#include <string>
#include <list>
#include "base_fmi.hpp"
#include "sdsl_fmi.hpp"
#include "fmi.hpp"
#include "timer.hpp"

void test(FMI &fmi, const std::vector<Base> &bases) {
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

    std::vector<Base> fwd_bases = seq_to_bases(fwd_str);

    BaseFMI base_fmi;
    SdslFMI sdsl_fmi;
    
    if (argc > 2) {
        index_fname = std::string(argv[2]);
        std::ifstream index_file(index_fname);
        base_fmi = BaseFMI(index_file, 128);
        sdsl_fmi.construct(fwd_str);
        index_file.close();
        //sdsl_fmi = SdslFMI(index_fname);
    } else {
        base_fmi = BaseFMI(fwd_str + "$", 128);
        sdsl_fmi.construct(fwd_str);
    }
    std::cout << "Built\n";
    
    Timer t;
    t.reset();

    test(base_fmi, fwd_bases);

    std::cout << t.lap() << "\n";

    test(sdsl_fmi, fwd_bases);
    
    std::cout << t.lap() << "\n";
}
