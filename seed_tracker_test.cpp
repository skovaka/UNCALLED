#include "seed_tracker.hpp"
#include <iostream>
#include <set>

bool check_aln(SeedTracker &t, std::string label, int min_len, int min_ratio) {
    if (t.longest_seed < min_len) {
        return false;
    }

    std::vector<ReadAln> alns = t.get_alignments(1);
    if (t.top_ratio(alns) < min_ratio) {
        return false;
    }

    ReadAln best = alns[0];

    std::cout << label << "\t" 
              << best.ref_st_.start_ << "-" 
              << best.ref_en_ << "\t"
              << best.evt_st_ << "-"
              << best.evt_en_ << "\t"
              << best.total_len_ << "\t" 
              << t.top_ratio(alns) << "\n";

    return true;
}

int main(int argc, char **argv) {

    if (argc < 3) {
        std::cerr << "Error: must include min length and ratio\n";
        return 1;
    }

    int min_len = atoi(argv[1]);
    double min_ratio = atof(argv[2]);
    
    std::string prev_strand, strand, evt_str, ref_str;
    double prob;

    int i, evt_st, evt_en, ref_st, ref_en;

    #ifdef DEBUG_PROB
    double min_evt_prob_;
    #endif

    SeedTracker fwd_tracker, rev_tracker;

    std::string time;

    bool read_aligned = false, parsing_read = false;

    std::string read_name;

    int seed_count = 0;

    while (std::cin) {

        std::cin >> strand;
        
        if (strand[0] == '=') {

            if (parsing_read) {
                parsing_read = false;

                if (!read_aligned) {
                    std::cout << read_name << " FAILED " << seed_count << "\n";
                }
                fwd_tracker.reset();
                rev_tracker.reset();
                
                read_aligned = false;
                seed_count = 0;
            } else {
                std::cin >> read_name;
                int i = read_name.find_last_of('/');
                if (i >= 0 && i < read_name.size()) {
                    read_name = read_name.substr(i+1);
                }
                parsing_read = true;
            }


            getline(std::cin, strand);
            continue;
        }

        std::cin >> time >> evt_str >> ref_str >> prob;

        seed_count++;

        #ifdef DEBUG_PROB
        std::cin >> min_evt_prob_;
        #endif

        i = evt_str.find('-');
        evt_st = atoi(evt_str.substr(0, i).c_str());
        evt_en = atoi(evt_str.substr(i + 1).c_str());

        i = ref_str.find('-');
        ref_st = atoi(ref_str.substr(0, i).c_str());
        ref_en = atoi(ref_str.substr(i + 1).c_str());

        Result r(evt_st, evt_en-evt_en, prob, ref_st, ref_en);

        if (strand == "fwd") {
            fwd_tracker.add_seed(r);

            if (!read_aligned && check_aln(fwd_tracker, read_name+" fwd "+time, min_len, min_ratio)) {
                read_aligned = true;
            }
        } else {
            rev_tracker.add_seed(r);

            if (!read_aligned && check_aln(rev_tracker, read_name+" rev "+time, min_len, min_ratio)) {
                read_aligned = true;
            }
        }

        #ifdef DEBUG_PROB
        std::cerr << ref_en << "\t" << min_evt_prob_ << "\n";
        #endif

    }
}
