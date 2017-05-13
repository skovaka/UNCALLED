#include <iostream>
#include <fstream>
#include <algorithm>
#include <tuple>
#include <vector>
#include <math.h>
#include <stddef.h>
#include "timer.h"
#include "nano_bwt.hpp"
#include "boost/math/distributions/students_t.hpp"


#define MER_LEN 6

const int alph_size = (int) pow(4, MER_LEN);
int BASE_IDS[(int) ('t'+1)];
mer_id mer_to_id(std::string seq, int i);
void init_base_ids();


/* TODO: read the actual file instead of the dummy version */
NanoFMI::NanoFMI(std::ifstream &model_in, std::vector<mer_id> &mer_seq, int tally_sp) {

    double em_mean, em_stdev, es_mean, es_stdev;
    
    std::cerr << "Reading model\n";

    for (int i = 0; i < alph_size; i++) {
        model_in >> em_mean >> em_stdev >> es_mean >> es_stdev;
        em_means.push_back(em_mean);
        em_stdevs.push_back(em_stdev);
        es_means.push_back(es_mean);
        es_stdevs.push_back(es_stdev);
    }

    
    tally_dist = tally_sp;
    
    mer_seq_tmp = &mer_seq;

    // create bwt
    make_bwt();
}

NanoFMI::NanoFMI(std::vector<double> model_in, std::vector<mer_id> &mer_seq, int tally_sp) {

    
    std::cerr << "Reading model\n";

    for (int i = 0; i < 4096; i++) {
        em_means.push_back(model_in[4*i+0]);
        em_stdevs.push_back(model_in[4*i+1]);
        es_means.push_back(model_in[4*i+2]);
        es_stdevs.push_back(model_in[4*i+3]);
    }

    em_means.push_back(-1);
    em_stdevs.push_back(-1);
    es_means.push_back(-1);
    es_stdevs.push_back(-1);

    tally_dist = tally_sp;
    
    mer_seq_tmp = &mer_seq;

    // create bwt
    make_bwt();
}

void NanoFMI::make_bwt() 
{

    //suffix_ar.resize(mer_seq_tmp->size());

    Timer timer;

    std::vector<unsigned int> suffix_ar_tmp(mer_seq_tmp->size());

    for (unsigned int i = 0; i < suffix_ar_tmp.size(); i++) {
        suffix_ar_tmp[i] = i;
    }


    std::cerr << "SA init time: " << timer.lap() << "\n";

    std::sort(suffix_ar_tmp.begin(), suffix_ar_tmp.end(), *this);

    suffix_ar.swap(suffix_ar_tmp);

    std::cerr << "SA sort time: " << timer.lap() << "\n";

    //Allocate space
    bwt.resize(mer_seq_tmp->size());
    mer_f_starts.resize(alph_size);
    mer_counts.resize(alph_size);
    mer_tally.resize(alph_size);

    for (mer_id i = 0; i < alph_size; i++)
        mer_tally[i].resize((mer_seq_tmp->size() / tally_dist) + 1);
    
    std::cerr << "FM init time: " << timer.lap() << "\n";

    int tally_mod = tally_dist;
    
    for (unsigned int i = 0; i < suffix_ar.size(); i++) {

        if (suffix_ar[i] > 0)
            bwt[i] = mer_seq_tmp->at(suffix_ar[i]-1);
        else
            bwt[i] = mer_seq_tmp->at(suffix_ar[suffix_ar.size()-1]);


        mer_counts[bwt[i]]++;

        if (tally_mod == tally_dist) {
            for (mer_id j = 0; j < alph_size; j++)
                mer_tally[j][i / tally_dist] = mer_counts[j];
            tally_mod = 0;
        }
        tally_mod += 1;
    }

    std::cerr << "FM build time: " << timer.lap() << "\n";

    for (mer_id i = 0; i < alph_size; i++) {
        mer_f_starts[i] = 1;
        for (mer_id j = 0; j < alph_size; j++)
            if (signal_compare(i, j) > 0)
                mer_f_starts[i] += mer_counts[j];
    }



    //for (mer_id i = 0; i < alph_size; i++)
    //    std::cout << i << "\t" << mer_counts[i] << "\n";
}
    
bool NanoFMI::operator() (unsigned int rot1, unsigned int rot2) {
    
    int c;
    for (unsigned int i = 0; i < mer_seq_tmp->size(); i++) {
        
        //if (rot1+i >= mer_seq_tmp->size())
        //    return true;
        //
        //if (rot2+i >= mer_seq_tmp->size())
        //    return false;

        c = signal_compare(mer_seq_tmp->at(rot1+i), 
                           mer_seq_tmp->at(rot2+i));
        
        if (c == 0)
            continue;

        if (c < 0)
            return true;

       return false;
    }

    return false;
}

int NanoFMI::signal_compare(mer_id mer1, mer_id mer2) {
    //if (mer1 >= em_means.size()) {
    //    if (mer2 >= em_means.size())
    //        return 0;
    //    else
    //        return -1;
    //} else if (mer2 >= em_means.size()) {
    //    return 1;
    //}

    if (em_means[mer1] < em_means[mer2])
        return -1;
    else if (em_means[mer1] > em_means[mer2])
        return 1;

    if (es_means[mer1] < es_means[mer2])
        return -1;
    else if (es_means[mer1] > es_means[mer2])
        return 1;

    if (em_stdevs[mer1] < em_stdevs[mer2])
        return 1;
    else if (em_stdevs[mer1] > em_stdevs[mer2])
        return -1;

    if (es_stdevs[mer1] < es_stdevs[mer2])
        return 1;
    else if (es_stdevs[mer1] > es_stdevs[mer2])
        return -1;

    return 0;
}

double NanoFMI::t_test(Event e, int i, ScaleParams scale) 
{
    // scale model mean
    double model_mean = em_means[i] * scale.scale + scale.shift;
    double model_stdv = em_stdevs[i] * scale.var;

    // calculate t-stat
    // https://en.wikipedia.org/wiki/Student%27s_t-test#Equal_or_unequal_sample_sizes.2C_equal_variance
    double t = (e.mean - model_mean) / (sqrt(2/e.length) * sqrt(pow(model_stdv, 2) + pow(e.stdv, 2)/2));
    int df = (e.length * 2) - 2;
    boost::math::students_t dist(df);
    double p = 2 * boost::math::cdf(dist, fabs(t));
    return p;
}

std::vector<MerRanges> NanoFMI::match_event(Event e, ScaleParams scale) {
    std::vector<MerRanges> mers;
    double stdv_scale = scale.var;
    for (mer_id i = 0; i < em_means.size(); i++) {
        if ((e.mean <= em_means[i] && e.mean + (e.stdv * stdv_scale) >= em_means[i] - (es_means[i] * stdv_scale)) ||
            (e.mean  > em_means[i] && e.mean - (e.stdv * stdv_scale) < em_means[i] + (es_means[i] * stdv_scale))) {
        //if (e.mean == em_means[i]) {
            MerRanges m = {i, std::vector<int>()};
            mers.push_back(m);
        }
    }
    return mers;
}

int NanoFMI::tally_cp_dist(int i) {
    int cp = (i / tally_dist)*tally_dist; //Closest checkpoint < i

    //Check if checkpoint after i is closer
    if (i - cp > (cp + tally_dist) - i && cp + tally_dist < bwt.size())
        cp += tally_dist;

    return cp - i;
}

int NanoFMI::get_tally(mer_id c, int i) {
    if (i < 0)
        return -1;

    int cp_dist = tally_cp_dist(i);
    int tally = mer_tally[c][(i + cp_dist) / tally_dist];


    if (cp_dist > 0) {
        for (unsigned int j = i+1; j <= i + cp_dist; j++)
            if (bwt[j] == c)
                tally--;

    } else if (cp_dist < 0) {
        for (unsigned int j = i; j > i + cp_dist; j--)
            if (bwt[j] == c)
                tally++;
    }

    return tally;
}

void NanoFMI::lf_map(std::vector<Event> events, ScaleParams scale) {
    std::vector<MerRanges> f_locs, l_locs;

    f_locs = match_event(events.back(), scale);
    for (unsigned int f = 0; f < f_locs.size(); f++) {
        f_locs[f].ranges.push_back(mer_f_starts[f_locs[f].mer]);
        f_locs[f].ranges.push_back(mer_f_starts[f_locs[f].mer]+mer_counts[f_locs[f].mer]-1);
    }

    bool mer_matched = false;

    for (int i = events.size()-2; i >= 0; i--) {

        l_locs = match_event(events[i], scale);
        mer_matched = false;

        int size = 0;
        for(unsigned int f = 0; f < f_locs.size(); f ++) 
            size += f_locs[f].ranges.size();

        //std::cout << "Aligning " << f_locs.size() << " " << size << " " << l_locs.size() << "\n";

        for (unsigned int l = 0; l < l_locs.size(); l++) { 
            for(unsigned int f = 0; f < f_locs.size(); f ++) {
                for (unsigned int r = 0; r < f_locs[f].ranges.size(); r += 2) {
                    int min_tally = get_tally(l_locs[l].mer, f_locs[f].ranges[r] - 1),
                        max_tally = get_tally(l_locs[l].mer, f_locs[f].ranges[r + 1]);
                    
                    //std::cout << f_locs[f].ranges[r] << " " << f_locs[f].ranges[r + 1] << "\n";
                    //std::cout << min_tally << " " << max_tally << "\n";

                    if (min_tally != max_tally) {
                        l_locs[l].ranges.push_back(mer_f_starts[l_locs[l].mer]+min_tally);
                        l_locs[l].ranges.push_back(mer_f_starts[l_locs[l].mer]+max_tally-1);
                        mer_matched = true;
                    }
                }
            }
        }

        f_locs.swap(l_locs);
        //l_locs.clear();

        if (!mer_matched)
            break;
    }

    if (mer_matched) {
        for (int l = 0; l < f_locs.size(); l++) {
            for (int r = 0; r < f_locs[l].ranges.size(); r += 2) {
                for (int s = f_locs[l].ranges[r]; s <= f_locs[l].ranges[r+1]; s++) {
                    std::cerr << "Match: " << suffix_ar[s] << "\n";
                }
            }
        }
    }

    std::cerr << mer_matched << "\n";
}



void init_base_ids() 
{
    for (unsigned int i = 0; i < ('t'+1); i++)
        BASE_IDS[i] = -1;
    BASE_IDS[(unsigned int) 'A'] = BASE_IDS[(unsigned int) 'a'] = 0;
    BASE_IDS[(unsigned int) 'C'] = BASE_IDS[(unsigned int) 'c'] = 1;
    BASE_IDS[(unsigned int) 'G'] = BASE_IDS[(unsigned int) 'g'] = 2;
    BASE_IDS[(unsigned int) 'T'] = BASE_IDS[(unsigned int) 't'] = 3;
}

//std::string id_to_mer(int id);
//def i_to_mer(i):
//    ret = [0]*MER_LEN
//    for n in range(0, MER_LEN):
//        ret[MER_LEN-n-1] = BASES[(i >> n*2) & 3]
//    return "".join(ret)

mer_id mer_to_id(std::string seq, unsigned int i) {
    mer_id id = BASE_IDS[(unsigned int) seq[i]];
    for (unsigned int j = 1; j < MER_LEN; j++)
        id = (id << 2) | BASE_IDS[(unsigned int) seq[i+j]];
    return id;
}

std::vector<mer_id> parse_fasta(std::ifstream &fasta_in) {
    init_base_ids();
    std::vector<mer_id> ids;
    std::string cur_line, prev_line, full_seq;

    getline(fasta_in, cur_line); //read past header

    while (getline(fasta_in, cur_line)) {
        if (prev_line.size() > 0)
            full_seq = prev_line.substr(prev_line.size()-MER_LEN+1) + cur_line;
        else
            full_seq = cur_line;
        
        for (unsigned int i = 0; i < full_seq.size() - MER_LEN + 1; i++)
            ids.push_back(mer_to_id(full_seq, i));

        prev_line = cur_line;
    }

    ids.push_back((int) pow(4, MER_LEN));

    return ids;
}



/*
int main(int argc, char **argv) {
    
    if (argc < 3) {
        std::cerr << "Usage: nano_bwt <pore_model.txt> <reference.fasta>\n";
        return 1;
    }

    init_base_ids();

    std::ifstream model_in(argv[1]), fasta_in(argv[2]);

    int suff_gap = atoi(argv[3]), tally_gap = atoi(argv[4]);

    std::vector<mer_id> ids = parse_fasta(fasta_in);

    NanoFMI bwt(model_in, ids, suff_gap, tally_gap);
}
*/
