#include <iostream>
#include <math.h>

#include "normalizer.hpp"
#include "pore_model.hpp"
#include "fast5_reader.hpp"
#include "bwa_index.hpp"
#include "dtw.hpp"

const std::string CONF_DIR(std::getenv("UNCALLED_CONF")),
                  DEF_MODEL = CONF_DIR + "/r94_5mers.txt",
                  DEF_CONF = CONF_DIR + "/defaults.toml";

//bool load_conf(int argc, char** argv, Conf &conf);
//

const KmerLen KLEN = KmerLen::k5;

int main(int argc, char** argv) {
    std::string index_prefix = std::string(argv[1]), 
                fast5_fname = std::string(argv[2]),
                chr = std::string(argv[3]);
    u32 st = atoi(argv[4]), en = atoi(argv[5]);
    std::string read_name = std::string(argv[6]);

    float diag_coef = 1,
          vert_coef = 1,
          horz_coef = 2;

    bool row_substr = true,
         col_substr = false,
         create_events = true;
         //aln_rev = false;


    DTWSubSeq substr = DTWSubSeq::NONE;
    if (row_substr && col_substr) {
        std::cerr << "Error: only use one substr\n";
        return 1;
    }

    if (row_substr) {
        substr = DTWSubSeq::ROW;
    } else if (col_substr) {
        substr = DTWSubSeq::COL;
    }

    const PoreModel<KLEN> model(DEF_MODEL, false);
    Normalizer norm(model.get_mean(), model.get_stdv());

    EventDetector::PRMS = event_detection_defaults;
    EventDetector evdt;

    BwaIndex<KLEN> idx(index_prefix);
    idx.load_pacseq();
    std::vector<u16> kmers = idx.get_kmers(chr, st, en);

    Fast5Reader fast5s_;
    fast5s_.add_fast5(fast5_fname);
    fast5s_.add_read(read_name);
    fast5s_.fill_buffer();
    ReadBuffer read = fast5s_.pop_read();

    auto costfn = [&model](float e, u16 k) {return -model.match_prob(e,k);};


    if (create_events) {
        norm.set_signal(evdt.get_means(read.get_raw()));
    } else {
        norm.set_signal(read.get_raw());
    }

    std::vector<float> signal;
    signal.reserve(norm.unread_size());
    while (!norm.empty()) signal.push_back(norm.pop());

    //Takes up too much space :(
    if (signal.size() > 500000) {
        return 1;
    }

    std::cerr << "Aligning\n";

    DTW<float, u16, decltype(costfn)>::Prms dtwp = {substr, diag_coef, vert_coef, horz_coef, costfn};

    DTW<float, u16, decltype(costfn)> fwd_dtw(signal, kmers, dtwp);

    //if (aln_rev) {
    //    std::cerr << "Aligning rev\n";
    //    DTW<float, u16, decltype(costfn)> rev_dtw(signal, rev_kmers, costfn, substr, 
    //                                              diag_coef, vert_coef, horz_coef);
    //    std::cout << fwd_dtw.mean_score() << "\t" << rev_dtw.mean_score() << "\n";
    //    if (fwd_dtw.mean_score() > rev_dtw.mean_score()) {
    //        //rev_dtw.print_path();
    //    }
    //        return 0;
    //}

    //fwd_dtw.print_path();
    std::cout << fwd_dtw.mean_score() << "\t"
              << signal.size() << "\t" 
              << kmers.size() << "\n";
    return 0;
}
