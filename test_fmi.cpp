#include <iostream>
#include <string>
#include <list>
#include <cstdlib>
#include "base_fmi.hpp"
#include "sdsl_fmi.hpp"
#include "fmi.hpp"
#include "timer.hpp"

void test_rank_select(FMI &fmi, const std::vector<Base> &bases) {
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

std::vector<size_t> test_sa(FMI &fmi, size_t num_locs) {
    std::vector<size_t> locs(num_locs);
    for (size_t i = 0; i < num_locs; i++) {
        locs[i] = fmi.sa(rand() % fmi.size());
    }
    return locs;
}

int main(int argc, char** argv) {
    //:std::string ref_fname(argv[1]), index_fname, fwd_str, rev_str;

    //std::ifstream ref_file(ref_fname);
    //parse_fasta(ref_file, fwd_str, rev_str, false);

    //std::vector<Base> fwd_bases = seq_to_bases(fwd_str);
    //
    //if (argc > 2) {
    //    //std::ifstream index_file(index_fname);
    //    //base_fmi = BaseFMI(index_file, 128);
    //    //sdsl_fmi.construct(fwd_str);
    //    //index_file.close();
    //} else {
    //    //base_fmi = BaseFMI(fwd_str + "$", 128);
    //    sdsl_fmi.construct(fwd_str);
    //}
    //std::cout << "Built\n";

    //BaseFMI base_fmi;
    std::string index_fname = std::string(argv[1]);
    SdslFMI sdsl_fmi(index_fname);

    //std::ifstream index_file(index_fname);
    //BaseFMI base_fmi(index_fname, 128);

    size_t num_locs = atoi(argv[2]);

    Timer t;
    t.reset();
    test_sa(sdsl_fmi, num_locs);
    std::cout << t.lap() << "\n";

    //test_rank_select(base_fmi, fwd_bases);
    //test_sa(base_fmi, locs);

    //std::cout << t.lap() << "\n";

    //test_rank_select(sdsl_fmi, fwd_bases);
    
}
