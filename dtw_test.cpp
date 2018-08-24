#include <iostream>
#include <math.h>

#include "kmer_model.hpp"
#include "fast5.hpp"
#include "arg_parse.hpp"
#include "timer.hpp"
#include "dtw.hpp"

double **full_dtw(const KmerModel &model,
                  const std::vector<Event> &read_events, 
                  const std::vector<u16> &ref_kmers,
                  bool local = true,
                  bool prob = true) {

    const double INF = 100000000;

    const int rd_len = read_events.size(),
              rf_len = ref_kmers.size();

    double **dtw_mat = new double * [rf_len];
    for (int i = 0; i < rf_len; i++) {
        dtw_mat[i] = new double[rd_len];
    }


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
                cost = -model.event_match_prob(read_events[j], ref_kmers[i]);
            } else {
                cost = fabs(model.lv_means_[ref_kmers[i]] - read_events[j].mean);
            }


            dtw_mat[i][j] = fmin(diag + cost, 
                            fmin(left + cost, 
                                 up + 2*cost));
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

bool get_raw(std::string filename, std::vector<float> &raw) {
    if (!fast5::File::is_valid_file(filename)) {
        std::cerr << "Error: '" << filename << "' is not a valid file \n";
    }

    try {
        fast5::File file;
        file.open(filename);
        
        if (!file.is_open()) {  
            std::cerr << "Error: unable to open '" 
                      << filename << "'\n";
            return false;
        }

        raw = file.get_raw_samples();
        return true;
        
    } catch (hdf5_tools::Exception& e) {
        std::cerr << "Error: hdf5 exception '" << e.what() << "'\n";
        return false;
    }

    return false;

}

bool get_events(std::string filename, std::vector<Event> &events, EventDetector &ed) {
    std::vector<float> raw;
    if (!get_raw(filename, raw)) {
        return false;
    }

    events = ed.get_all_events(raw);
    return true;
}
    
enum Opt {MODEL     = 'm',
          REF_FASTA = 'x',
          READ_LIST = 'r'};

int main(int argc, char** argv) {

    ArgParse args("DTW tester");
    args.add_string(Opt::MODEL, "model", "/home-4/skovaka1@jhu.edu/code/nanopore_aligner/kmer_models/r9.4_180mv_450bps_6mer/5mers_pre.txt", "Nanopore kmer model");
    args.add_string(Opt::REF_FASTA, "ref_fasta", "", "");
    args.add_string(Opt::READ_LIST, "read_pattern", "", "");
    args.parse_args(argc, argv);

    const KmerModel model(args.get_string(Opt::MODEL), false);

    detector_param params = event_detection_defaults;
    params.threshold1 = 1.4;
    params.threshold2 = 1.1;
    EventDetector event_detector(params, 30, 150);

    std::ifstream ref_in(args.get_string(Opt::REF_FASTA));
    std::vector<u16> fwd_kmers, rev_kmers;
    model.parse_fasta(ref_in, fwd_kmers, rev_kmers);

    std::ifstream reads_in(args.get_string(Opt::READ_LIST));

    std::string read_fname;

    while (getline(reads_in, read_fname)) {

        //std::vector<Event> events;
        std::vector<float> events;
        //if (!get_events(read_fname, events, event_detector)) {
        if (!get_raw(read_fname, events)) {
            std::cout << "#" << read_fname << "\n"
                      << "-1\t-1\n#END#\n";
            continue;
        }

        model.normalize_raw(events);

        //bool full_read = false;

        auto cost = [&model](float e, u16 k) {return -model.event_match_prob({e,1,4},k);};
        DTW<float, u16, decltype(cost)> 
            fwd_dtw(events, fwd_kmers, cost, SubSeqDTW::ROW, 20, 1000000, 1),
            rev_dtw(events, rev_kmers, cost, SubSeqDTW::ROW, 20, 1000000, 1);

        std::cout << "== fwd: " << fwd_dtw.mean_score()
                  << " rev: " << rev_dtw.mean_score() << " ==\n";

        auto &ref_dtw = fwd_dtw.mean_score() < rev_dtw.mean_score()
                                   ? fwd_dtw : rev_dtw;

        ref_dtw.print_path();
    }

        /*
        double **fwd_mat = raw_dtw(model, 
                                   fwd_kmers, raw, 
                                   true, full_read),

               **rev_mat = raw_dtw(model, 
                                   rev_kmers, raw,
                                   true, full_read);
        
        double **best_mat = fwd_mat;
        bool fwd_strand = fwd_mat[fwd_kmers.size()-1][raw.size()-1]
                          < rev_mat[rev_kmers.size()-1][raw.size()-1];
        if (!fwd_strand) {
            best_mat = rev_mat;
        }

        std::vector<u16> &ref_kmers = fwd_strand ? fwd_kmers : rev_kmers;

        std::list< std::pair<size_t, size_t> > path 
            = get_path(best_mat, 
                       ref_kmers.size(), 
                       raw.size(),
                       true,
                       full_read);

        std::cout << "== " << read_fname << " ==\n";
        std::cout << "== " << (fwd_strand ? "fwd" : "rev") << " ==\n";

        double prev_prob = 0;
        for (auto p = path.begin(); p != path.end(); p++) {
            std::cout << p->second << "\t" //read coord
                      << p->first << "\t"  //ref coord
                      << raw[p->second] << "\t"
                      << ref_kmers[p->first] << "\t"
                      << best_mat[p->first][p->second] << "\n";
                      
            //std::cout << p->second << "\t" //read coord
            //          << p->first << "\t"  //ref coord
            //          << read_events[p->second].mean << "\t"
            //          << read_events[p->second].stdv << "\t"
            //          << read_events[p->second].length << "\t"
            //          << ref_kmers[p->first] << "\t"
            //          << prev_prob - best_mat[p->first][p->second] << "\n";
            
            prev_prob = best_mat[p->first][p->second];
        }

        std::cout << "== END ==\n";
    */
}
