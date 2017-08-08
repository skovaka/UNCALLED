/* sigalign.cpp 
 * align oxford nanopore event-level signals directly to the reference without
 * basecalling */

#include <iostream>
#include <string>
#include <unordered_map>
#include <math.h>

#include "fast5.hpp"
#include "kmer_model.hpp"
#include "nano_fmi.hpp"
#include "timer.h"


void align_kmers(std::string name, std::string strand, NanoFMI& fmi, std::vector<Event>& events, NormParams norm)
{
    Timer timer;
    int k = 16;

    // align each k-mer of length k
    //for (int i = events.size()-1; i >= k; i--) {
    for (int i = 2266; i >= 2226; i--) { //Should map to 12198-12214 - 55 events -> 36 genome BP
        auto matches = fmi.lf_map(events, i, k, norm);

        for (auto match = matches.begin(); match != matches.end(); match++) {
            std::cout << strand << "\t" 
                      << match->qry_start << "-" << match->qry_end << "\t"
                      << match->ref_start << "-" << match->ref_end << "\t"
                      << match->prob << "\n";
        }

        break;
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

    std::cerr << "Loading model\n";
    KmerModel model(model_fname);

    std::cerr << "Parsing fasta\n";
    std::vector<mer_id> fwd_ids, rev_ids;
    std::ifstream ref_file(ref_fname);
    model.parse_fasta(ref_file, fwd_ids, rev_ids);

    Timer timer;

    std::cerr << "Building forward FMI\n";
    //NanoFMI fwd_fmi(model, fwd_ids, tally_gap);

    std::cerr << "Building reverse FMI\n";
    NanoFMI rev_fmi(model, rev_ids, tally_gap);

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

            if (f.have_eventdetection_events()) {
                events = f.get_eventdetection_events();
            } else {
                std::cerr << "file " << argv[fix] << " does not contain events. skipping..." << std::endl;
                continue;
            }

            NormParams scale = model.get_norm_params(events);

            align_kmers(argv[fix], "rev", rev_fmi, events, scale);
            //align_kmers(argv[fix], "fwd", fwd_fmi, events, scale);

        }
        catch (hdf5_tools::Exception& e)
        {
            std::cerr << "hdf5 error: " << e.what() << std::endl;
        }
    }
}
