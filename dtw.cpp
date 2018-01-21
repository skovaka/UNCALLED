#include <iostream>
#include <math.h>

#include "kmer_model.hpp"
#include "fast5.hpp"
#include "arg_parse.hpp"
#include "timer.h"

double **full_dtw(const KmerModel &model,
                  const std::vector<Event> &read_events, 
                  const std::vector<mer_id> &ref_kmers,
                  const NormParams &norm,
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

bool get_events(std::string filename, std::vector<Event> &events) {
    events.clear();

    if (!fast5::File::is_valid_file(filename)) {
        std::cerr << "Error: '" << filename << "' is not a valid file \n";
        return false;
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
    
enum Opt {MODEL        = 'm',
          REF_FASTA    = 'i',
          READ_PATTERN = 'r',
          LOCAL        = 'l',
          PROB         = 'p'};

int main(int argc, char** argv) {

    ArgParse args("DTW tester");

    args.add_string(Opt::MODEL, "model", "/home-4/skovaka1@jhu.edu/code/nanopore_aligner/kmer_models/r9.4_180mv_450bps_6mer/5mers_pre.txt", "Nanopore kmer model");
    args.add_string(Opt::REF_FASTA,   "ref_fasta", "", "");
    args.add_string(Opt::READ_PATTERN,   "read_pattern", "", "");
    //args.add_string(Opt::REF, "reference", "", "Reference to align to (fasta)");
    args.add_flag(Opt::LOCAL, "local", "Will perform local alignment if provided, otherwise global");
    args.add_flag(Opt::PROB, "prob", "Will use kmer/event match probability for cost if probability, otherwise absolute mean difference");

    args.parse_args(argc, argv);

    std::cout << "#" << args.get_string(Opt::MODEL) << "\n"
              << "#" << args.get_string(Opt::REF_FASTA) << "\n";

    KmerModel model(args.get_string(Opt::MODEL));

    std::string read_prefix = args.get_string(Opt::READ_PATTERN), read_suffix;
    int i = read_prefix.find("*");
    read_suffix = read_prefix.substr(i+1);
    read_prefix = read_prefix.substr(0, i);

    std::ifstream ref_file(args.get_string(Opt::REF_FASTA));

    if (!ref_file) {
        std::cerr << "Error: couldn't open '" 
                  << args.get_string(Opt::REF_FASTA) << "'\n";
        return 1;
    }

    std::string ref_header, ref_seq, read_filename;

    while (getline(ref_file, ref_header)) {

        getline(ref_file, ref_seq);

        read_filename = read_prefix + ref_header.substr(1) + read_suffix;

        std::vector<Event> read_events;
        if (!get_events(read_filename, read_events)) {
            std::cout << "#" << read_filename << "\n"
                      << "-1\t-1\n#END#\n";
            continue;
        }

        NormParams norm = model.get_norm_params(read_events);

        std::vector<mer_id> ref_kmers(ref_seq.size() - model.kmer_len() + 1);
        for (size_t i = 0; i < ref_kmers.size(); i++) {
            ref_kmers[i] = model.kmer_to_id(ref_seq, i);
        }

        double **dtw_mat = full_dtw(model, 
                                    read_events, 
                                    ref_kmers, 
                                    norm,
                                    args.get_flag(Opt::LOCAL),
                                    args.get_flag(Opt::PROB));

        std::list< std::pair<size_t, size_t> > path 
            = get_path(dtw_mat, 
                       ref_kmers.size(), 
                       read_events.size(),
                       args.get_flag(Opt::LOCAL));
        
        std::cout << "#" << read_filename << "\n"
                  << "#" << norm.scale << " " << norm.shift << "\n";

        double prev_prob = 0;
        for (auto p = path.begin(); p != path.end(); p++) {
            std::cout << p->second << "\t" 
                      << p->first << "\t"
                      << read_events[p->second].mean << "\t"
                      << read_events[p->second].stdv << "\t"
                      << read_events[p->second].length << "\t"
                      << ref_kmers[p->first] << "\t"
                      << dtw_mat[p->first][p->second] - prev_prob << "\n";
            
            prev_prob = dtw_mat[p->first][p->second];
        }

        std::cout << "#END#\n";


        for (size_t i = 0; i < ref_kmers.size(); i++) {
            delete[] dtw_mat[i];
        }
        delete[] dtw_mat;

       
    }
}
