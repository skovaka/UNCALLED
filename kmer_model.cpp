#include <fstream>
#include <sstream>
#include <numeric>
#include <tuple>
#include <string>
#include <cmath>
#include <stddef.h>
#include "boost/math/distributions/students_t.hpp"
#include "kmer_model.hpp"

#define PI 3.1415926535897

KmerModel::KmerModel(std::string model_fname) {
    std::ifstream model_in(model_fname);

    //Read header and count number of columns
    std::string header;
    std::getline(model_in, header);
    int num_cols = 0;
    bool prev_ws = true;
    for (unsigned int i = 0; i < header.size(); i++) {
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
    mer_id k_id;
    double lv_mean, lv_stdv, sd_mean, sd_stdv, lambda, weight;

    //Check if table includes "ig_lambda" column
    bool has_lambda = num_cols >= 7;

    //Read first line
    if (has_lambda) {
        model_in >> kmer >> lv_mean >> lv_stdv >> sd_mean 
                 >> sd_stdv >> lambda >> weight;
        lambda_ = lambda;
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
    sd_means_ = new double[kmer_count_+1];
    sd_stdvs_ = new double[kmer_count_+1];

    lv_means_[kmer_count_] = 
    lv_vars_x2_[kmer_count_] = 
    lognorm_denoms_[kmer_count_] = 
    sd_means_[kmer_count_] = 
    sd_stdvs_[kmer_count_] = -1;

    model_mean_ = 0;

    //Stores which kmers can follow which
    neighbors_ = new std::list<mer_id>[kmer_count_]; 
    rev_comp_ids_ = new mer_id[kmer_count_]; 

    //Read and store rest of the model
    do {

        //Get unique ID for the kmer
        k_id = kmer_to_id(kmer);

        //Check if kmer is valid
        if (k_id < 0 || k_id >= kmer_count_) {
            std::cerr << "Error: kmer '" << kmer << "' is invalid\n";

        //Store kmer information
        } else {

            //Store kmer reverse complement
            rev_comp_ids_[k_id] = kmer_to_id(reverse_complement(kmer));

            //Store all neighboring kmers
            for (short b = 0; b < 4; b++) { //4 magic number?
                neighbors_[k_id].push_back(
                    kmer_to_id(id_to_base(b) + kmer.substr(0, k_ - 1)));
            }

            //Store model information
            lv_means_[k_id] = lv_mean;
            lv_vars_x2_[k_id] = 2*lv_stdv*lv_stdv;
            lognorm_denoms_[k_id] = log(sqrt(PI * lv_vars_x2_[k_id]));

            sd_means_[k_id] = sd_mean;
            sd_stdvs_[k_id] = sd_stdv;

            model_mean_ += lv_mean;
        }

        //Read next line
        if (has_lambda) {
            model_in >> kmer >> lv_mean >> lv_stdv >> sd_mean 
                     >> sd_stdv >> lambda >> weight;
            
            //Lambda changes - might as well use sd_stdv
            if (lambda_ != lambda)
                lambda_ = -1;
        } else {
            model_in >> kmer >> lv_mean >> lv_stdv 
                     >> sd_mean >> sd_stdv >> weight;
        }

    //Read until eof or correct number of kmers read
    } while (!model_in.eof());
    
    //Compute model level mean and stdv
    model_mean_ /= kmer_count_;
    model_stdv_ = 0;
    for (mer_id k_id = 0; k_id < kmer_count_; k_id++)
        model_stdv_ += pow(lv_means_[k_id] - model_mean_, 2);
    model_stdv_ = sqrt(model_stdv_ / kmer_count_);
}


KmerModel::~KmerModel() {
    delete [] lv_means_;
    delete [] lv_vars_x2_;
    delete [] lognorm_denoms_;
    delete [] sd_means_;
    delete [] sd_stdvs_;
    delete [] neighbors_; 
    delete [] rev_comp_ids_; 
}

bool KmerModel::event_valid(const Event &e) const {
    return e.mean > 0 && e.stdv >= 0 && e.length > 0;
}

float KmerModel::event_match_prob(Event e, mer_id k_id, NormParams norm) const {
    //double norm_mean = lv_means_[k_id] * norm.scale + norm.shift;
    //return (-pow(e.mean - norm_mean, 2) / lv_vars_x2_[k_id]) - lognorm_denoms_[k_id];

    double norm_mean = (e.mean - norm.shift) / norm.scale;
    return (-pow(norm_mean - lv_means_[k_id], 2) / lv_vars_x2_[k_id]) - lognorm_denoms_[k_id];
}

float KmerModel::get_stay_prob(Event e1, Event e2) const {

    if (e1.length == 0 || e2.length == 0)
        return 0;

    if (e1 == e2)
        return log(1);

    double var1 = e1.stdv*e1.stdv, var2 = e2.stdv*e2.stdv;

    double t = (e1.mean - e2.mean) / sqrt(var1/e1.length + var2/e2.length);

    int df = pow(var1/e1.length + var2/e2.length, 2) 
               / (pow(var1/e1.length, 2) / (e1.length-1) 
                    + pow(var2/e2.length, 2) / (e2.length-1));

    boost::math::students_t dist(df);
    double q = boost::math::cdf(boost::math::complement(dist, fabs(t)));

    return log(q);
}

std::string KmerModel::reverse_complement(std::string &seq) const {
    std::string rev(seq.size(), 'N');
    for (unsigned int i = 0; i < seq.size(); i++) {
        char c = 'N';
        switch(seq[i]) {
            case 'A':
            case 'a':
            c = 'T';
            break;
            case 'T':
            case 't':
            c = 'A';
            break;
            case 'G':
            case 'g':
            c = 'C';
            break;
            case 'C':
            case 'c':
            c = 'G';
            break;
        }
        rev[rev.size()-i-1] = c;
    }

    return rev;
}

mer_id KmerModel::kmer_to_id(std::string kmer, int offset) const {
    mer_id id = base_to_id(kmer[offset]);
    for (unsigned int j = 1; j < k_; j++)
        id = (id << 2) | base_to_id(kmer[offset+j]);
    return id;
}
    
short KmerModel::base_to_id(char b) {
    switch (b) {
        case 'A':
        case 'a':
        return 0;

        case 'C':
        case 'c':
        return 1;

        case 'G':
        case 'g':
        return 2;

        case 'T':
        case 't':
        return 3;
    }

    return -1;
}

char KmerModel::id_to_base(short i) {
    switch (i) {
        case 0:
        return 'A';
        
        case 1:
        return 'C';

        case 2:
        return 'G';

        case 3:
        return 'T';
    }

    return 'N';
}

NormParams KmerModel::get_norm_params(const std::vector<Event> &events) const {

    //Compute events mean
    double events_mean = 0;
    for (auto e = events.begin(); e != events.end(); e++) 
        events_mean += e->mean;
    events_mean /= events.size();

    //Compute events stdv
    double events_stdv = 0;
    for (auto e = events.begin(); e != events.end(); e++) 
        events_stdv += pow(e->mean - events_mean, 2);
    events_stdv = sqrt(events_stdv / events.size());
    
    /* get scaling parameters */
    NormParams params;
    params.scale = events_stdv / model_stdv_;
    params.shift = events_mean - (params.scale * model_mean_);

    return params;
}

//Reads the given fasta file and stores forward and reverse k-mers
void KmerModel::parse_fasta(
                 std::ifstream &fasta_in, 
                 std::vector<mer_id> &fwd_ids, 
                 std::vector<mer_id> &rev_ids) const {

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
        for (unsigned int i = 0; i < fwd_seq.size() - k_ + 1; i++)
            fwd_ids.push_back(kmer_to_id(fwd_seq, i));
    
        prev_line = cur_line;
    }

    rev_ids = std::vector<mer_id>(fwd_ids.size());
    for (unsigned int i = 0; i < fwd_ids.size(); i++)
        rev_ids[rev_ids.size()-i-1] = rev_comp_ids_[fwd_ids[i]];

    //Add EOF character
    fwd_ids.push_back((int) kmer_count_);
    rev_ids.push_back((int) kmer_count_);

}
