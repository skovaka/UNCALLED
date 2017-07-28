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


void align_kmers(std::string name, std::string strand, NanoFMI& fmi, std::vector<Event>& events, ScaleParams scale)
{
    Timer timer;
    int k = 16;
    int aligned_kmers = 0;

    // align each k-mer of length k
    //for (int i = events.size()-1; i >= k; i--) {
    for (int i = 14831; i >= k; i--) {
        //std::vector<Event> kmer(events.begin() + i, events.begin() + i + k);
        //int i = 14830;
        int count  = 0;
        //std::cout << "Aligning " << i << "\n";
        count = fmi.lf_map(events, i, k, scale);
        if (count) {
            std::cout << strand << "\t" <<  k << "\t" << timer.lap() << "\t" << i << "\t" << count << std::endl;
            //aligned_kmers += 1;
        }
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

    /* navigate through all the read files provided */
    for (int fix = 4; fix < argc; fix++) {
        std::vector<Event> events;
        if (not fast5::File::is_valid_file(argv[fix])) {
            std::cerr << "<" << argv[fix] << "> is not a valid file. skipping... " << std::endl;
            continue;
        }
        fast5::File f;
        try
        {
            f.open(argv[fix]);
            assert(f.is_open());
            // std::cerr << "processing read " << argv[fix] << std::endl;
            //make sure that the fast5 contains an event sequence 
            if (f.have_eventdetection_events()) {
                events = f.get_eventdetection_events();
            } else {
                std::cerr << "file " << argv[fix] << " does not contain events. skipping..." << std::endl;
                continue;
            }
            ScaleParams scale = get_scale_params(model, events);
            //TODO: start alignment, pass in scale, model (?), events
            align_kmers(argv[fix], "rev", rev_fmi, events, scale);
            align_kmers(argv[fix], "fwd", fwd_fmi, events, scale);
            //int fwd_count = fwd_fmi.lf_map(events, scale);
            //int rev_count = rev_fmi.lf_map(events, scale);
            //std::cout << argv[fix] << ": " << fwd_count << " matches on fwd strand, " << rev_count << " matches on rev strand" << std::endl;
        }
        catch (hdf5_tools::Exception& e)
        {
            std::cerr << "hdf5 error: " << e.what() << std::endl;
        }
    }
}
