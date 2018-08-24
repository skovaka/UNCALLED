#ifndef _INCL_MODEL_TOOLS
#define _INCL_MODEL_TOOLS

#include <array>
#include <utility>
#include <list>
#include "basepairs.hpp"
#include "event_detector.hpp"

#define NUM_BASES 4

typedef unsigned short int u16;

typedef struct NormParams {
    double shift, scale;
} NormParams;

class KmerModel {
    public:
    typedef std::list<u16>::const_iterator neighbor_itr;

    KmerModel(std::string model_fname, bool complement);
    ~KmerModel();

    NormParams get_norm_params(const std::vector<Event> &events) const;
    void normalize_events(std::vector<Event> &events, NormParams norm={0, 0}) const;
    void normalize_raw(std::vector<float> &raw) const;

    std::pair<neighbor_itr, neighbor_itr> get_neighbors(u16 k_id) const {
        return std::pair<neighbor_itr, neighbor_itr>
                    (neighbors_[k_id].begin(), neighbors_[k_id].end());
    }

    bool event_valid(const Event &e) const;

    float event_match_prob(Event evt, u16 k_id) const;

    float get_stay_prob(Event e1, Event e2) const; 

    inline size_t kmer_len() const {return k_;}
    inline size_t kmer_count() const {return kmer_count_;}

    u16 kmer_to_id(std::string kmer, int offset = 0) const;
    u16 kmer_comp(u16 kmer);
    //std::string id_to_kmer(u16 kmer) const;

    u8 get_first_base(u16 k) const;
    u8 get_last_base(u16 k) const;
    u8 get_base(u16 kmer, size_t i) const;

    void parse_fasta (std::ifstream &fasta_in, 
                 std::vector<u16> &fwd_ids, 
                 std::vector<u16> &rev_ids) const;

    double *lv_means_, *lv_vars_x2_, *lognorm_denoms_, *sd_means_, *sd_stdvs_;

    private:
    unsigned short k_, kmer_count_;
    double lambda_, model_mean_, model_stdv_;
    bool complement_;

    std::list<u16> *neighbors_; //TODO: NO MORE LISTS
    u16 *rev_comp_ids_;
};

#endif

