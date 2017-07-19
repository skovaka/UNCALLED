#ifndef INCL_NANOBWT
#define INCL_NANOBWT

#include "model_tools.hpp"
#include <list>

void parse_fasta(std::ifstream &fasta_in, 
                 std::vector<mer_id> &fwd_ids, 
                 std::vector<mer_id> &rev_ids);

//TODO: make nested class?
//class MerRanges {
//    public:
//
//    MerRanges(){};
//    MerRanges(int start, int range_len, float prob);
//    MerRanges(MerRanges prev, float prob);
//
//    void add_match(int start, int range_len, float prob, int match_len);
//    void add_stay(MerRanges prev, float prob);
//
//    std::vector<int> starts, ends, match_lens;
//    std::vector<float> prob_sums;
//};

class FMQuery {
    public:
    static const char NEXT = 'n', STAY = 't', SKIP = 'k';
    int start, end, match_len, stays;
    float prob_sum;
    std::list<char> align;

    FMQuery(){}

    FMQuery(int range_start, int range_len, float prob)
        : start(range_start), end(range_start+range_len-1), 
          match_len(1), stays(0), prob_sum(prob) {
        align.push_front(NEXT);          
    }

    FMQuery(FMQuery prev, int range_start, int range_len, float prob)
        : start(range_start), end(range_start+range_len-1), 
          match_len(prev.match_len+1), stays(prev.stays), prob_sum(prev.prob_sum+prob), align(prev.align) {
        align.push_front(NEXT);          
    }

    FMQuery(FMQuery prev, float prob)
        : start(prev.start), end(prev.end), 
          match_len(prev.match_len), stays(prev.stays+1), prob_sum(prev.prob_sum+prob), align(prev.align) {
        align.push_front(STAY);
    }
};

class NanoFMI 
{
    public:

    NanoFMI(std::ifstream &model_in, std::vector<mer_id> &mer_seq, int tally_sp);

    NanoFMI(std::vector<double> model_in, std::vector<mer_id> &mer_seq, int tally_sp);
    //std::vector<mer_id> get_bwt();

    bool operator() (unsigned int rot1, unsigned int rot2);

    int lf_map(std::vector<Event> events, int seed_end, int match_len, ScaleParams scale);

    private:
    int signal_compare(mer_id mer1, mer_id mer2);
    int tally_cp_dist(int i);
    int get_tally(mer_id c, int i);
    float get_evt_prob(Event e, int i, ScaleParams scale);
    float get_stay_prob(Event e1, Event e2);
    std::vector<double> em_means, em_stdevs, es_means, es_stdevs;
    std::vector<mer_id> *mer_seq_tmp;
    
    //Sparseness of suffix and tally arrays
    int tally_dist;

    std::vector<mer_id> bwt;                        //L array
    std::vector<int> mer_counts, mer_f_starts;      //F array
    std::vector<unsigned int> suffix_ar;            //Suffix array
    std::vector< std::vector<int> > mer_tally;      //Rank tally


    void make_bwt();
    
};

#endif

