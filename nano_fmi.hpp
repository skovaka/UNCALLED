#ifndef INCL_NANOBWT
#define INCL_NANOBWT

#include "kmer_model.hpp"
#include <list>

void parse_fasta(std::ifstream &fasta_in, 
                 std::vector<mer_id> &fwd_ids, 
                 std::vector<mer_id> &rev_ids);

class FMQuery {
    public:
    //const char NEXT = 'n', STAY = 't', SKIP = 'k';
    int start, end, match_len, stays;
    float prob_sum;
    //std::list<char> align;

    FMQuery(){ std::cout << "NOPE\n";}

    FMQuery(const FMQuery &fq){
        start = fq.start;
        end = fq.end;
        match_len = fq.match_len;
        stays = fq.stays; 
        prob_sum = fq.prob_sum; 
    }

    FMQuery(int range_start, int range_len, float prob)
        : start(range_start), end(range_start+range_len-1), 
          match_len(1), stays(0), prob_sum(prob) {
    }

    FMQuery(const FMQuery &prev, int range_start, int range_len, float prob) {
        start = range_start;
        end = range_start+range_len-1;
        match_len = prev.match_len+1;
        stays = prev.stays; 
        prob_sum = prev.prob_sum+prob; 
    }

    FMQuery(const FMQuery &prev, float prob)
        : start(prev.start), end(prev.end), 
          match_len(prev.match_len), stays(prev.stays+1), prob_sum(prev.prob_sum+prob) {
        
    }

    bool same_range(const FMQuery &q) const;

    float avg_prob() const;

    friend bool operator< (const FMQuery &q1, const FMQuery &q2);
};

class NanoFMI 
{
    public:
    NanoFMI(KmerModel &model, std::vector<mer_id> &mer_seq, int tally_dist);

    bool operator() (unsigned int rot1, unsigned int rot2);

    int lf_map(std::vector<Event> &events, int seed_end, 
               int match_len, NormParams norm_params);

    private:
    int signal_compare(mer_id mer1, mer_id mer2);
    int tally_cp_dist(int i);
    int get_tally(mer_id c, int i);

    float get_stay_prob(Event e1, Event e2); //hmm

    KmerModel *model_;
    std::vector<mer_id> *mer_seq_;
    
    //Sparseness of suffix and tally arrays
    int tally_dist_;

    std::vector<mer_id> bwt_;                        //L array
    std::vector<int> mer_counts_, mer_f_starts_;      //F array
    std::vector<unsigned int> suffix_ar_;            //Suffix array
    std::vector< std::vector<int> > mer_tally_;      //Rank tally


    void make_bwt();
    
};

#endif

