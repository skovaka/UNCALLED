#ifndef _INCL_MODEL_TOOLS
#define _INCL_MODEL_TOOLS

#include <array>
#include <utility>
#include <list>
#include "basepairs.hpp"
#include "fast5.hpp"

#define NUM_BASES 4

typedef unsigned char base_t;
typedef unsigned short int mer_id;
typedef fast5::EventDetection_Event Event;
typedef fast5::Basecall_Event BCEvent;

typedef struct NormParams {
    double shift, scale;
} NormParams;

class KmerModel {
    public:
    typedef std::list<mer_id>::const_iterator neighbor_itr;

    KmerModel(std::string model_fname);
    ~KmerModel();

    NormParams get_norm_params(const std::vector<Event> &events) const;
    void normalize_events(std::vector<Event> &events, NormParams norm={0, 0}) const;

    std::pair<neighbor_itr, neighbor_itr> get_neighbors(mer_id k_id) const {
        return std::pair<neighbor_itr, neighbor_itr>
                    (neighbors_[k_id].begin(), neighbors_[k_id].end());
    }

    bool event_valid(const Event &e) const;

    float event_match_prob(Event evt, mer_id k_id) const;

    float get_stay_prob(Event e1, Event e2) const; 

    inline size_t kmer_len() const {return k_;}
    inline size_t kmer_count() const {return kmer_count_;}

    mer_id kmer_to_id(std::string kmer, int offset = 0) const;
    //std::string id_to_kmer(mer_id kmer) const;

    base_t get_first_base(mer_id k) const;
    base_t get_last_base(mer_id k) const;
    base_t get_base(mer_id kmer, size_t i) const;

    void parse_fasta (std::ifstream &fasta_in, 
                 std::vector<mer_id> &fwd_ids, 
                 std::vector<mer_id> &rev_ids) const;

    double *lv_means_, *lv_vars_x2_, *lognorm_denoms_, *sd_means_, *sd_stdvs_;

    private:
    unsigned short k_, kmer_count_;
    double lambda_, model_mean_, model_stdv_;

    std::list<mer_id> *neighbors_; 
    mer_id *rev_comp_ids_;
};

#endif

