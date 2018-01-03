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

    NormParams norm = {0, 1};//= model.get_norm_params(read_events);

    double cost, up, left, diag;

    for (int i = 0; i < rf_len; i++) {
        for (int j = 0; j < rd_len; j++) {

            up   = i > 0 ? dtw_mat[i-1][j] : INF;

            left = j > 0 ? dtw_mat[i][j-1] : (local ? 0 :  INF);

            diag = j > 0 && i > 0 ? dtw_mat[i-1][j-1] : (j == i || (j == 0 && local) ? 0 : INF);

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
    
int main(int argc, char** argv) {
    bool local = true;

    KmerModel model(argv[1]);

    std::vector<Event> read_events;
    read_events.push_back( {85.083612, 1.517846, 5} );
    read_events.push_back( {85.083612, 1.517846, 5} );
    read_events.push_back( {76.635809, 1.705015, 5} );
    read_events.push_back( {76.635809, 1.705015, 5} );
    read_events.push_back( {76.635809, 1.705015, 5} );
    read_events.push_back( {87.429392, 1.962027, 5} );
    read_events.push_back( {106.44446, 2.191321, 5} );
    read_events.push_back( {106.44446, 2.191321, 5} );
    read_events.push_back( {106.44446, 2.191321, 5} );
    read_events.push_back( {91.063935, 2.709875, 5} );
    read_events.push_back( {91.063935, 2.709875, 5} );
    read_events.push_back( {85.299201, 1.517846, 5} );
    read_events.push_back( {85.083612, 1.517846, 5} );
    read_events.push_back( {76.635809, 1.705015, 5} );
    read_events.push_back( {87.429392, 1.962027, 5} );
    read_events.push_back( {106.44446, 2.191321, 5} );
    read_events.push_back( {91.063935, 2.709875, 5} );
    read_events.push_back( {85.299201, 1.517846, 5} );

    std::vector<mer_id> ref_kmers;
    ref_kmers.push_back(0);
    ref_kmers.push_back(1);
    ref_kmers.push_back(4);
    ref_kmers.push_back(16);
    ref_kmers.push_back(64);
    ref_kmers.push_back(256);
    ref_kmers.push_back(0);
    ref_kmers.push_back(0);
    ref_kmers.push_back(0);
    ref_kmers.push_back(0);
    ref_kmers.push_back(1);
    ref_kmers.push_back(4);
    ref_kmers.push_back(16);
    ref_kmers.push_back(64);
    ref_kmers.push_back(256);

    double **dtw_mat = full_dtw(model, read_events, ref_kmers, local, true);

    for (size_t i = 0; i < ref_kmers.size(); i++) {
        for (size_t j = 0; j < read_events.size(); j++) {
            std::cout << dtw_mat[i][j] << "\t";
        }
        std::cout << "\n";
    }

    std::list< std::pair<size_t, size_t> > path = get_path(dtw_mat, 
                                                           ref_kmers.size(), 
                                                           read_events.size(),
                                                           local);
    
    for (auto p = path.begin(); p != path.end(); p++) {
        std::cout << p->first << "\t" << p->second << "\n";
    }

}
