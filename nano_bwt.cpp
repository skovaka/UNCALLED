#include <iostream>
#include <fstream>
#include <algorithm>
#include <unordered_map>
#include <vector>
#include <math.h>
#include "timer.h"
#include "nano_bwt.hpp"

#define MER_LEN 6

const int alph_size = (int) pow(4, MER_LEN);
int BASE_IDS[(int) ('t'+1)];
mer_id mer_to_id(std::string seq, int i);
void init_base_ids();


/* TODO: read the actual file instead of the dummy version */
NanoFMI::NanoFMI(std::ifstream &model_in, std::vector<mer_id> mer_seq_tmp, int sa_sp, int tally_sp) {

    double em_mean, em_stdev, es_mean, es_stdev;
    
    std::cerr << "Reading model\n";

    for (int i = 0; i < alph_size; i++) {
        model_in >> em_mean >> em_stdev >> es_mean >> es_stdev;
        em_means.push_back(em_mean);
        em_stdevs.push_back(em_stdev);
        es_means.push_back(es_mean);
        es_stdevs.push_back(es_stdev);
    }

    sa_dist = sa_sp;
    tally_dist = tally_sp;
    
    // mer_seq_tmp = &mer_seq;
    mer_seq = mer_seq_tmp;

    // create bwt
    make_bwt();
}

NanoFMI::NanoFMI(std::vector<double> model_in, std::vector<mer_id> mer_seq_tmp, int sa_sp, int tally_sp) {

    
    std::cerr << "Reading model\n";

    for (int i = 0; i < 4096; i++) {
        em_means.push_back(model_in[4*i+0]);
        em_stdevs.push_back(model_in[4*i+1]);
        es_means.push_back(model_in[4*i+2]);
        es_stdevs.push_back(model_in[4*i+3]);
    }

    sa_dist = sa_sp;
    tally_dist = tally_sp;
    
    //mer_seq_tmp = &mer_seq;
    mer_seq = mer_seq_tmp;
    // create bwt
    make_bwt();
}

void NanoFMI::make_bwt() 
{

    std::vector<unsigned int> suffix_ar(mer_seq.size());

    Timer timer;

    for (unsigned int i = 0; i < suffix_ar.size(); i++) {
        suffix_ar[i] = i;
    }

    std::cerr << "SA init time: " << timer.lap() << "\n";

    std::sort(suffix_ar.begin(), suffix_ar.end(), *this);

    std::cerr << "SA sort time: " << timer.lap() << "\n";

    //Allocate space
    bwt.resize(mer_seq.size());
    mer_counts.resize(alph_size);
    mer_tally.resize(alph_size);
    sparse_sa.reserve((mer_seq.size() / sa_dist) + 1);
    
    for (mer_id i = 0; i < alph_size; i++)
        mer_tally[i].resize((mer_seq.size() / tally_dist) + 1);
    
    std::cerr << "FM init time: " << timer.lap() << "\n";

    int tally_mod = 0;
    
    for (unsigned int i = 0; i < suffix_ar.size(); i++) {
        if (suffix_ar[i] % sa_dist == 0)
            sparse_sa[i] = suffix_ar[i];

        if (suffix_ar[i] > 0)
            bwt[i] = mer_seq[suffix_ar[i]-1];
        else
            bwt[i] = mer_seq[suffix_ar[suffix_ar.size()-1]];

        mer_counts[bwt[i]]++;

        if (tally_mod == tally_dist) {
            for (mer_id j = 0; j < alph_size; j++)
                mer_tally[j][i / tally_dist] = mer_counts[bwt[i]];
            tally_mod = 0;
        }
        tally_mod += 1;
    }

    std::cerr << "FM build time: " << timer.lap() << "\n";

    //for (mer_id i = 0; i < alph_size; i++)
    //    std::cout << i << "\t" << mer_counts[i] << "\n";
}
    
bool NanoFMI::operator() (unsigned int rot1, unsigned int rot2) {
    
    int c;
    for (unsigned int i = 0; i < mer_seq.size(); i++) {
        
        if (rot1+i >= mer_seq.size())
            return true;
        
        if (rot2+i >= mer_seq.size())
            return false;

        c = signal_compare(mer_seq.at(rot1+i), 
                           mer_seq.at(rot2+i));
        
        if (c == 0)
            continue;

        if (c < 0)
            return true;

       return false;
    }

    return false;
}

int NanoFMI::signal_compare(mer_id mer1, mer_id mer2) {
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

    //ids.push_back((int) pow(4, MER_LEN));

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
