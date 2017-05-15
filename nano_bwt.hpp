#ifndef INCL_NANOBWT
#define INCL_NANOBWT

#include "model_tools.hpp"

void parse_fasta(std::ifstream &fasta_in, 
                 std::vector<mer_id> &fwd_ids, 
                 std::vector<mer_id> &rev_ids);


typedef struct MerRanges {
    mer_id mer;
    std::vector<int> ranges;
} MerRanges;

class NanoFMI 
{
    public:

    NanoFMI(std::ifstream &model_in, std::vector<mer_id> &mer_seq, int tally_sp);

    NanoFMI(std::vector<double> model_in, std::vector<mer_id> &mer_seq, int tally_sp);
    //std::vector<mer_id> get_bwt();

    bool operator() (unsigned int rot1, unsigned int rot2);

    void lf_map(std::vector<Event> events, ScaleParams scale);

    private:
    int signal_compare(mer_id mer1, mer_id mer2);
    std::vector<MerRanges> match_event(Event e, ScaleParams scale);
    int tally_cp_dist(int i);
    int get_tally(mer_id c, int i);
    double t_test(Event e, int i, ScaleParams scale);
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

