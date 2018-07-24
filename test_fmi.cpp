#include <iostream>
#include <string>
#include <list>
#include <cstdlib>
#include "base_fmi.hpp"
#include "bwa_fmi.hpp"
#include "sdsl_fmi.hpp"
#include "bwa/bwt.h"
#include "fmi.hpp"
#include "timer.hpp"

void test_rank_select(FMI &fmi, const std::vector<u8> &bases) {
    Range r = fmi.get_full_range(bases.back());

    for (size_t i = bases.size()-2; i < bases.size(); i--) {
        r = fmi.get_neighbor(r, bases[i]);
    }

    for (size_t s = r.start_; s <= r.end_; s++) {
        std::cout << s << " " << fmi.sa(s) << "\n";
    }
}

std::vector<size_t> test_sa(FMI &fmi, size_t num_locs) {
    std::vector<size_t> locs(num_locs);
    for (size_t i = 0; i < num_locs; i++) {
        locs[i] = fmi.sa(rand() % fmi.size());
    }
    return locs;
}

//Probably not really useful
//Just how the BWT is stored, not queried
std::vector<u32> str_to_bytes(std::string sseq) {
    std::vector<u32> bseq(sseq.size() / 16 + (sseq.size() % 16 != 0));
    u32 i = 0;
    for (u32 j = 0; j < bseq.size(); j++) {
        for (u8 k = 30; k < 32 && i < sseq.size(); k -= 2) {
            bseq[j] |= (BASE_BYTES[(u32) sseq[i++]] << k);
        }
    }
    return bseq;
}

//Probably not really useful
//Just how the BWT is stored, not queried
std::string bytes_to_str(std::vector<u32> bseq, u32 len) {
    std::string sseq(len, 'A');
    u32 i = 0;
    for (u32 j = 0; j < bseq.size(); j++) {
        for (u8 k = 30; k < 32 && i < sseq.size(); k -= 2) {
            sseq[i++] = BASE_CHARS[(bseq[j] >> k) & 3];
        }
    }
    return sseq;
}

int main(int argc, char** argv) {
    //BaseFMI base_fmi;
    std::string index_prefix = std::string(argv[1]);

    std::cout << "Converting query...\n";
    std::string query_fname(argv[2]);
    std::ifstream query_in(query_fname);
    std::string fwd_squery, rev_squery;
    parse_fasta(query_in, fwd_squery, rev_squery, false);
    //std::string squery = std::string(argv[2]);
    std::vector<u8> bquery(fwd_squery.size());
    for (u32 i = 0; i < fwd_squery.size(); i++) {
        bquery[i] = BASE_BYTES[(u8) fwd_squery[i]];
    }
    std::cout << "Query length: " << fwd_squery.size() << "\n\n";

    std::string bwt_fname = index_prefix + ".bwt",
                sa_fname = index_prefix + ".sa",
                base_fname = index_prefix + ".ufmi";

    Timer t;

    std::cout << "Loading BWA...\n";
    t.reset();
    BwaFMI bwa_fmi(index_prefix);
    std::cout << "Time: " << t.get() << ", Size: " << bwa_fmi.size() << "\n\n";

    std::cout << "Loading UNCALLED...\n";
    t.reset();
    std::ifstream in(index_prefix + ".ufmi");
    BaseFMI base_fmi(in, 128);
    std::cout << "Time: " << t.get() << ", Size: " << base_fmi.size() << "\n\n";

    std::cout << "Loading SDSL...\n";
    t.reset();
    SdslFMI sdsl_fmi(index_prefix + ".sfmi");
    std::cout << "Time: " << t.get() << ", Size: " << sdsl_fmi.size() << "\n\n";

    size_t NUM_SA = 1000000;

    t.reset();
    test_sa(bwa_fmi, NUM_SA);
    std::cout << "BWA SA time: " << t.get() << "\n\n";

    t.reset();
    test_sa(base_fmi, NUM_SA);
    std::cout << "UNCALLED SA time: " << t.get() << "\n\n";
    
    t.reset();
    test_sa(sdsl_fmi, NUM_SA);
    std::cout << "SDSL SA time: " << t.get() << "\n\n";

    t.reset();
    test_rank_select(bwa_fmi, bquery);
    std::cout << "BWA R/S time: " << t.get() << "\n\n";

    t.reset();
    test_rank_select(base_fmi, bquery);
    std::cout << "UNCALLED R/S time: " << t.get() << "\n\n";

    t.reset();
    test_rank_select(sdsl_fmi, bquery);
    std::cout << "SDSL R/S time: " << t.get() << "\n\n";
}
