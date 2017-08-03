#ifndef _INCL_MODEL_TOOLS
#define _INCL_MODEL_TOOLS

#include <array>
#include "fast5.hpp"

#define NUM_BASES 4

typedef unsigned short int mer_id;
typedef fast5::EventDetection_Event Event;

typedef struct NormParams {
    double shift, scale;
} NormParams;

class KmerModel {
    public:
    //KmerModel(std::ifstream model_in);

    KmerModel(std::string model_fname);

    NormParams get_norm_params(std::vector<Event> events);

    std::array<mer_id, NUM_BASES>::const_iterator neighbor_begin(mer_id k_id) {
        return neighbors_[k_id].begin();
    }

    std::array<mer_id, NUM_BASES>::const_iterator neighbor_end(mer_id k_id) {
        return neighbors_[k_id].end();
    }

    float event_match_prob(Event evt, mer_id k_id, NormParams norm);

    int compare_kmers(mer_id k_id1, mer_id k_id2);

    inline int kmer_len() {return k_;}
    inline int kmer_count() {return kmer_count_;}

    void parse_fasta(std::ifstream &fasta_in, 
                 std::vector<mer_id> &fwd_ids, 
                 std::vector<mer_id> &rev_ids);
    private:
    unsigned short k_, kmer_count_;
    double *lv_means_, *lv_stdvs_, *sd_means_, *sd_stdvs_;
    double lambda_, model_mean_, model_stdv_;

    std::array<mer_id, NUM_BASES> *neighbors_; 
    mer_id *rev_comp_ids_;


    mer_id kmer_to_id(std::string kmer, int offset = 0);
    std::string reverse_complement(std::string &seq);

    static char id_to_base(short i);
    static short base_to_id(char b);
};

#endif

