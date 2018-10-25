#include <fstream>
#include <sstream>
#include <numeric>
#include <tuple>
#include <string>
#include <cmath>
#include <stddef.h>
//#include "boost/math/distributions/students_t.hpp"
#include "kmer_model.hpp"

#define PI 3.1415926535897

std::string get_reverse_complement(const std::string &seq) {
    std::string rev(seq.size(), 'N');
    for (u64 i = 0; i < seq.size(); i++) {
        rev[rev.size()-i-1] = BASE_COMP_C[(u8) seq[i]];
    }
    return rev;
}

std::string get_complement(const std::string &seq) {
    std::string comp(seq.size(), 'N');
    for (u64 i = 0; i < seq.size(); i++) {
        comp[i] = BASE_COMP_C[(u8) seq[i]];
    }
    return comp;
}

KmerModel::KmerModel(std::string model_fname, bool complement) {
    std::ifstream model_in(model_fname);

    //Read header and count number of columns
    std::string header;
    std::getline(model_in, header);
    u8 num_cols = 0;
    bool prev_ws = true;
    for (u8 i = 0; i < header.size(); i++) {
        if (header[i] == ' ' || header[i] == '\t') {
            prev_ws = true;
        } else {
            if (prev_ws)
                num_cols++;
            prev_ws = false;
        }
    }

    //Variables for reading model
    std::string kmer, neighbor_kmer;
    u16 k_id;
    double lv_mean, lv_stdv, sd_mean, sd_stdv, lambda, weight;

    //Check if table includes "ig_lambda" column
    bool has_lambda = num_cols >= 7;

    //Read first line
    if (has_lambda) {
        model_in >> kmer >> lv_mean >> lv_stdv >> sd_mean 
                 >> sd_stdv >> lambda >> weight;
        lambda_ = lambda;
    } else if (num_cols == 4) {
        model_in >> kmer >> lv_mean >> lv_stdv >> sd_mean; 
        lambda_ = -1;
    } else {
        model_in >> kmer >> lv_mean >> lv_stdv >> sd_mean 
                 >> sd_stdv >> weight;
        lambda_ = -1;
    }

    //Compute number of kmers (4^k) and reserve space for model
    k_ = kmer.size();
    kmer_count_ = pow(4, k_); //4 magic number?
    lv_means_ = new double[kmer_count_+1];
    lv_vars_x2_ = new double[kmer_count_+1];
    lognorm_denoms_ = new double[kmer_count_+1];

    lv_means_[kmer_count_] = 
    lv_vars_x2_[kmer_count_] = 
    lognorm_denoms_[kmer_count_] = -1;

    model_mean_ = 0;

    //Stores which kmers can follow which
    neighbors_.resize(kmer_count_);
    rev_comp_ids_ = new u16[kmer_count_]; 

    //Read and store rest of the model
    do {
        //Complement kmer if needed
        if (complement) {
            kmer = get_complement(kmer);
        }

        //Get unique ID for the kmer
        k_id = kmer_to_id(kmer);


        //Check if kmer is valid
        if (k_id < 0 || k_id >= kmer_count_) {
            std::cerr << "Error: kmer '" << kmer << "' is invalid\n";

        //Store kmer information
        } else {

            //Store kmer reverse complement
            rev_comp_ids_[k_id] = kmer_to_id(get_reverse_complement(kmer));

            //Store all neighboring kmers
            for (short b = 0; b < 4; b++) { //4 magic number?
                //neighbors_[k_id].push_back(kmer_to_id(BASE_CHARS[b] + kmer.substr(0, k_ - 1)));
                neighbors_[k_id].push_back(kmer_to_id(kmer.substr(1, k_) + BASE_CHARS[b]));
           }

            //Store model information
            lv_means_[k_id] = lv_mean;
            lv_vars_x2_[k_id] = 2*lv_stdv*lv_stdv;
            lognorm_denoms_[k_id] = log(sqrt(PI * lv_vars_x2_[k_id]));

            model_mean_ += lv_mean;
        }

        //Read first line
        if (has_lambda) {
            model_in >> kmer >> lv_mean >> lv_stdv >> sd_mean 
                     >> sd_stdv >> lambda >> weight;
            lambda_ = lambda;
        } else if (num_cols == 4) {
            model_in >> kmer >> lv_mean >> lv_stdv >> sd_mean; 
            lambda_ = -1;
        } else {
            model_in >> kmer >> lv_mean >> lv_stdv >> sd_mean 
                     >> sd_stdv >> weight;
            lambda_ = -1;
        }

    //Read until eof or correct number of kmers read
    } while (!model_in.eof());
    
    //Compute model level mean and stdv
    model_mean_ /= kmer_count_;
    model_stdv_ = 0;
    for (u16 k_id = 0; k_id < kmer_count_; k_id++)
        model_stdv_ += pow(lv_means_[k_id] - model_mean_, 2);
    model_stdv_ = sqrt(model_stdv_ / kmer_count_);
}


KmerModel::~KmerModel() {
    delete [] lv_means_;
    delete [] lv_vars_x2_;
    delete [] lognorm_denoms_;
    delete [] rev_comp_ids_; 
}

bool KmerModel::event_valid(const Event &e) const {
    return e.mean > 0 && e.stdv >= 0 && e.length > 0;
}

float KmerModel::event_match_prob(float e, u16 k_id) const {
    return (-pow(e - lv_means_[k_id], 2) / lv_vars_x2_[k_id]) - lognorm_denoms_[k_id];
}

float KmerModel::event_match_prob(const Event &e, u16 k_id) const {
    return event_match_prob(e.mean, k_id);
}

u16 KmerModel::kmer_to_id(std::string kmer, u64 offset) const {
    u16 id = BASE_BYTES[(u8) kmer[offset]];
    for (u8 j = 1; j < k_; j++) {
        id = (id << 2) | BASE_BYTES[(u8) kmer[offset+j]];
    }
    return id;
}

u16 KmerModel::kmer_comp(u16 kmer) {
    return kmer ^ ((1 << (2*k_)) - 1);
}

u8 KmerModel::get_base(u16 kmer, u8 i) const {
    return (u8) ((kmer >> (2 * (k_-i-1))) & 0x3);
}

u8 KmerModel::get_first_base(u16 kmer) const {
    return (u8) ((kmer >> (2*k_ - 2)) & 0x3);
}

u8 KmerModel::get_last_base(u16 kmer) const {
    return (u8) (kmer & 0x3);
}

NormParams KmerModel::get_norm_params(const std::vector<Event> &events) const {
    //Compute events mean
    double events_mean = 0;
    for (auto e : events) 
        events_mean += e.mean;
    events_mean /= events.size();

    //Compute events stdv
    double events_stdv = 0;
    for (auto e : events) 
        events_stdv += pow(e.mean - events_mean, 2);
    events_stdv = sqrt(events_stdv / events.size());
    
    /* get scaling parameters */
    NormParams params;
    params.scale = model_stdv_ / events_stdv;
    params.shift = model_mean_ - (params.scale * events_mean);

    return params;
}

NormParams KmerModel::get_norm_params(const std::vector<float> &events) const {
    //Compute events mean
    double events_mean = 0;
    for (auto e : events) 
        events_mean += e;
    events_mean /= events.size();

    //Compute events stdv
    double events_stdv = 0;
    for (auto e : events) 
        events_stdv += pow(e - events_mean, 2);
    events_stdv = sqrt(events_stdv / events.size());
    
    /* get scaling parameters */
    NormParams params;
    params.scale = model_stdv_ / events_stdv;
    params.shift = model_mean_ - (params.scale * events_mean);

    return params;
}

void KmerModel::normalize(std::vector<Event> &events, NormParams norm) const {
    if (norm.scale == 0) {
        norm = get_norm_params(events);
    }
    for (size_t i = 0; i < events.size(); i++) {
        events[i].mean = norm.scale * events[i].mean + norm.shift;
    }
}

void KmerModel::normalize(std::vector<float> &events, NormParams norm) const {
    if (norm.scale == 0) {
        norm = get_norm_params(events);
    }
    for (size_t i = 0; i < events.size(); i++) {
        events[i] = norm.scale * events[i] + norm.shift;
    }
}

u16 KmerModel::get_neighbor(u16 k, u8 i) const {
    return neighbors_[k][i];
}

//Reads the given fasta file and stores forward and reverse k-mers
void KmerModel::parse_fasta(
                 std::ifstream &fasta_in, 
                 std::vector<u16> &fwd_ids, 
                 std::vector<u16> &rev_ids) const {

    //For parsing the file
    std::string cur_line, prev_line, fwd_seq;
    
    getline(fasta_in, cur_line); //read past header

    while (getline(fasta_in, cur_line)) {

        //Saves end of previous line to get last five 6-mers
        if (prev_line.size() > 0)
            fwd_seq = prev_line.substr(prev_line.size() - k_ + 1) + cur_line;
        else
            fwd_seq = cur_line;
        
        //Store forward IDs
        for (u64 i = 0; i < fwd_seq.size() - k_ + 1; i++)
            fwd_ids.push_back(kmer_to_id(fwd_seq, i));
    
        prev_line = cur_line;
    }

    rev_ids = std::vector<u16>(fwd_ids.size());
    for (u64 i = 0; i < fwd_ids.size(); i++)
        rev_ids[rev_ids.size()-i-1] = rev_comp_ids_[fwd_ids[i]];
}
