/* sigalign.cpp 
 * align oxford nanopore event-level signals directly to the reference without
 * basecalling */

#include <iostream>
#include <string>
#include <unordered_map>
#include <math.h>

#include "fast5.hpp"
#include "model_tools.hpp"
#include "nano_bwt.hpp"
#include "timer.h"


void align_kmers(NanoFMI& fmi, std::vector<Event>& events, ScaleParams scale)
{
    Timer timer;
    // for each k-mer length (step:k=2*k)
    for (int k = 8; k < pow(2, 5) + 1; k = 2*k) {
        int aligned_kmers = 0;
        // align each k-mer of length k
        for (int i = 0; i < events.size() - k + 1; i++) {
        // for (int i = 0; i < 1000; i++) {
            std::vector<Event> kmer(events.begin() + i, events.begin() + i + k);
            int count  = 0;
            count = fmi.lf_map(kmer, scale);
            if (count) {
                aligned_kmers += 1;
                std::cout << k << "\t" << timer.lap() << "\t" << i << "\t" << count << std::endl;
            }
        }
        // std::cout << "for " << events.size() - k + 1 << " " << k << "-mers: " << aligned_kmers << " aligned." << std::endl;
    }

}

int main(int argc, char** argv) 
{
    if (argc < 5) { 
        std::cerr << "need a file. " << 
            "(usage: sigalign reference model tally_sp fast5_file_1 [fast5_file_2, ...]" << std::endl;
        return 1;
    }
    std::string ref_fname = argv[1];
    std::string model_fname = argv[2];
    int tally_gap = std::stoi(argv[3]);

    /* load the model into a vector:
     * model has size 4096*4. every group of four represents a 6mer, and is
     * ordered lexocographically from AAAAAA to TTTTTT. format:
     * [level_mean1, level_stdv1, sd_mean1, sd_stdv1, 
     *  level_mean2, level_stdv2, sd_mean2, sd_stdv2, 
     *  ...]
     */
    std::vector<double> model = load_model(model_fname);

    /* initialize the FM-index */
    std::ifstream ref_fp(ref_fname);
    std::vector<mer_id> fwd_ids, rev_ids;

    Timer timer;

    std::cerr << "Parsing fasta\n";
    parse_fasta(ref_fp, fwd_ids, rev_ids);

    std::cerr << "Building forward FMI\n";
    NanoFMI fwd_fmi(model, fwd_ids, tally_gap);

    std::cerr << "Building reverse FMI\n";
    NanoFMI  rev_fmi(model, rev_ids, tally_gap);
    std::cout << "read_name\tk\ttime\tpos_in_read\tmatches" << std::endl;

    std::vector<Event> read = simulate_read(model, rev_ids, 0, 500);


    std::cerr << "Mapping to reverse index\n";
    align_kmers(rev_fmi, read, ScaleParams());

    std::cerr << "Mapping to foward index\n";
    align_kmers(fwd_fmi, read, ScaleParams());

    std::cerr << "Done\n";
}
