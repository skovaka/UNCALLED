#include <iostream>
#include <math.h>

#include "kmer_model.hpp"
#include "fast5.hpp"
#include "arg_parse.hpp"
#include "timer.h"

double **full_dtw(const KmerModel &model,
                  const std::vector<Event> &read_events, 
                  const std::vector<mer_id> &ref_kmers,
                  bool local = true,
                  bool prob = true) {

    const double INF = 100000000;

    const int rd_len = read_events.size(),
              rf_len = ref_kmers.size();

    double **dtw_mat = new double * [rf_len];
    for (int i = 0; i < rf_len; i++) {
        dtw_mat[i] = new double[rd_len];
    }

    NormParams norm = model.get_norm_params(read_events);

    double cost, up, left, diag;

    for (int i = 0; i < rf_len; i++) {
        for (int j = 0; j < rd_len; j++) {

            up   = i > 0 ? dtw_mat[i-1][j] 
                         : INF;

            left = j > 0 ? dtw_mat[i][j-1] 
                         : (local ? 0 
                                  : INF);

            diag = j > 0 && i > 0 ? dtw_mat[i-1][j-1] 
                                  : (j == i || (j == 0 && local) ? 0 
                                                                 : INF);

            if (prob) {
                cost = -model.event_match_prob(read_events[j], ref_kmers[i], norm);
            } else {
                cost = fabs(model.lv_means_[ref_kmers[i]] - read_events[j].mean);
            }

            dtw_mat[i][j] = fmin(diag, fmin(left, up)) + cost;
        }
    }

    return dtw_mat;
}

std::list< std::pair<size_t, size_t> > get_path(double **dtw_mat, 
                                                int rf_len, 
                                                int rd_len,
                                                bool local = true) {

    std::pair<size_t, size_t> p = {rf_len-1, rd_len-1};

    if (local) {
        for (int i = 0; i < rf_len; i++) {
            if (dtw_mat[i][rd_len-1] < dtw_mat[p.first][rd_len-1]) {
                p.first = i;
            } 
        }
    }

    std::list< std::pair<size_t, size_t> > path;
    path.push_front(p);

    double up, left, diag;

    bool done = false;
    while(!done) {
        if (p.first == 0) {
            p = {p.first, p.second-1};
        } else if (p.second == 0) {
            p = {p.first-1, p.second};
        } else {
            up = dtw_mat[p.first-1][p.second];
            left = dtw_mat[p.first][p.second-1];
            diag = dtw_mat[p.first-1][p.second-1];

            if (diag <= left && diag <= up) {
                p = {p.first-1, p.second-1};
            } else if (left <= up) {
                p = {p.first, p.second-1};
            } else {
                p = {p.first-1, p.second};
            }
        }

        path.push_front(p);
        
        done = p.second == 0 && (local || p.first == 0); 
    }

    return path;
}

bool get_bcevents(std::string filename, std::vector<BCEvent> &events) {
    events.clear();

    if (!fast5::File::is_valid_file(filename)) {
        std::cerr << "Error: '" << filename << "' is not a valid file \n";
    }

    try {
        fast5::File file;
        file.open(filename);
        
        if (!file.is_open()) {  
            std::cerr << "Error: unable to open '" 
                      << filename << "'\n";

        } else if (!file.have_basecall_events(0)) {
            std::cerr << "Error: file '" << filename
                      << "' does not contain events\n";

        } else {
            events = file.get_basecall_events(0);
            return true;
        }

    } catch (hdf5_tools::Exception& e) {
        std::cerr << "Error: hdf5 exception '" << e.what() << "'\n";
    }

    return false;

}

bool get_events(std::string filename, std::vector<Event> &events) {
    events.clear();

    if (!fast5::File::is_valid_file(filename)) {
        std::cerr << "Error: '" << filename << "' is not a valid file \n";
    }

    try {
        fast5::File file;
        file.open(filename);
        
        if (!file.is_open()) {  
            std::cerr << "Error: unable to open '" 
                      << filename << "'\n";

        } else if (!file.have_eventdetection_events()) {
            std::cerr << "Error: file '" << filename
                      << "' does not contain events\n";

        } else {
            events = file.get_eventdetection_events();
            return true;
        }

    } catch (hdf5_tools::Exception& e) {
        std::cerr << "Error: hdf5 exception '" << e.what() << "'\n";
    }

    return false;
}
    
enum Opt {MODEL  = 'm',
          READ   = 'i',
          REF    = 'x',
          LOCAL  = 'l',
          PROB   = 'p',
          REV    = 'r'};

int main(int argc, char** argv) {

    ArgParse args("DTW tester");

    args.add_string(Opt::MODEL, "model", "/home-4/skovaka1@jhu.edu/code/nanopore_aligner/kmer_models/r9.4_180mv_450bps_6mer/5mers_pre.txt", "Nanopore kmer model");
    args.add_string(Opt::READ, "read", "", "Event detected read (fast5)");
    args.add_string(Opt::REF, "reference", "", "Reference to align to (fasta)");
    args.add_flag(Opt::LOCAL, "local", "Will perform local alignment if provided, otherwise global");
    args.add_flag(Opt::PROB, "prob", "Will use kmer/event match probability for cost if probability, otherwise absolute mean difference");
    args.add_flag(Opt::REV, "rev", "Will align to reverse complement of reference if provided");

    args.parse_args(argc, argv);

    KmerModel model(args.get_string(Opt::MODEL));

    std::vector<mer_id> fwd_kmers, rev_kmers, ref_kmers;
    std::ifstream ref_file(args.get_string(Opt::REF));
    model.parse_fasta(ref_file, fwd_kmers, rev_kmers);

    if (args.get_flag(Opt::REV)) {
        ref_kmers = rev_kmers;
    } else  {
        ref_kmers = fwd_kmers;
    }

    ref_kmers.pop_back();

    std::vector<Event> read_events;
    if (!get_events(args.get_string(Opt::READ), read_events)) {
        return 1;
    }

    //std::vector<BCEvent> read_bcevents;
    //if (!get_bcevents(args.get_string(Opt::READ), read_bcevents)) {
    //    return 1;
    //}

    double **dtw_mat = full_dtw(model, 
                                read_events, 
                                ref_kmers, 
                                args.get_flag(Opt::LOCAL),
                                args.get_flag(Opt::PROB));

    //for (size_t i = ref_kmers.size()-11; i < ref_kmers.size(); i++) {
    //    for (size_t j = read_events.size()-11; j < read_events.size(); j++) {
    //        std::cout << dtw_mat[i][j] << "\t";
    //    }
    //    std::cout << "\n";
    //}

    //for (size_t i = 0; i < ref_kmers.size(); i++) {
    //    for (size_t j = 0; j < read_events.size(); j++) {
    //        std::cout << dtw_mat[i][j] << "\t";
    //    }
    //    std::cout << "\n";
    //}

    std::list< std::pair<size_t, size_t> > path = get_path(dtw_mat, 
                                                           ref_kmers.size(), 
                                                           read_events.size(),
                                                           args.get_flag(Opt::LOCAL));
    
    for (auto p = path.begin(); p != path.end(); p++) {
        std::cout << p->first << "\t" << p->second << "\t" << dtw_mat[p->first][p->second] << "\n";
    }

}
