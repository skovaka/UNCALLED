#ifndef INCL_NANOBWT
#define INCL_NANOBWT


typedef unsigned short int mer_id;
std::vector<mer_id> parse_fasta(std::ifstream &fasta_in);

class NanoFMI 
{
    public:

    NanoFMI(std::ifstream &model_in, std::vector<mer_id> mer_seq_tmp, int sa_sp, int tally_sp);

    NanoFMI(std::vector<double> model_in, std::vector<mer_id> mer_seq_tmp, int sa_sp, int tally_sp);
    //std::vector<mer_id> get_bwt();
    

    bool operator() (unsigned int rot1, unsigned int rot2);

    private:
    int signal_compare(mer_id mer1, mer_id mer2);

    std::vector<double> em_means, em_stdevs, es_means, es_stdevs;
    std::vector<mer_id> mer_seq;
    
    //Sparseness of suffix and tally arrays
    int tally_dist, sa_dist;

    std::vector<mer_id> bwt;                        //L array
    std::vector<int> mer_counts;                    //F array
    std::vector< std::vector<int> > mer_tally;      //Rank tally
    std::unordered_map<unsigned int, unsigned int> sparse_sa; //Sparse suffix array

    private:

    void make_bwt();
    
};

#endif

