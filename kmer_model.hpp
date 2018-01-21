#ifndef _INCL_MODEL_TOOLS
#define _INCL_MODEL_TOOLS

#include <array>
#include <utility>
#include <list>
#include "fast5.hpp"

#define NUM_BASES 4

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

    std::pair<neighbor_itr, neighbor_itr> get_neighbors(mer_id k_id) const {
        return std::pair<neighbor_itr, neighbor_itr>
                    (neighbors_[k_id].begin(), neighbors_[k_id].end());
    }

    bool event_valid(const Event &e) const;

    float event_match_prob(Event evt, mer_id k_id, NormParams norm) const;

    float get_stay_prob(Event e1, Event e2) const; 

    inline int kmer_len() const {return k_;}
    inline int kmer_count() const {return kmer_count_;}

    mer_id kmer_to_id(std::string kmer, int offset = 0) const;

    void parse_fasta (std::ifstream &fasta_in, 
                 std::vector<mer_id> &fwd_ids, 
                 std::vector<mer_id> &rev_ids) const;

    double *lv_means_, *lv_vars_x2_, *lognorm_denoms_, *sd_means_, *sd_stdvs_;

    private:
    unsigned short k_, kmer_count_;
    double lambda_, model_mean_, model_stdv_;

    std::list<mer_id> *neighbors_; 
    mer_id *rev_comp_ids_;


    std::string reverse_complement(std::string &seq) const ;

    static char id_to_base(short i);
    static short base_to_id(char b);
};

#endif

