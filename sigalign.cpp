/* sigalign.cpp 
 * align oxford nanopore event-level signals directly to the reference without
 * basecalling */

#include <iostream>
#include <string>
#include <unordered_map>
#include <math.h>

#include "fast5.hpp"
#include "seed_graph.hpp"
#include "kmer_model.hpp"
#include "nano_fmi.hpp"
#include "timer.h"


void align_kmers(const std::string &strand, 
                 KmerModel &model, 
                 NanoFMI& fmi, 
                 std::vector<Event>& events, 
                 NormParams norm) {
    int k = 32;

    std::cout << "== " << events.size() << " events ==\n";

    SeedGraph sg(model, fmi, norm, k, events.size() - k + 1, -5.29, -3.35, -8.343, 0.7);
    //SeedGraph sg(model, fmi, norm, k, 2266 - 2226 + 1, -9.2103, -3.75, -5.298, 0.7);

    // align each k-mer of length k
    for (int i = events.size()-1; i >= k; i--) {
    //for (int i = 2266; i >= 2226; i--) { //Should map to 12198-12214 - 55 events -> 36 genome BP
        std::vector<Result> results = sg.add_event(events[i]);

        for (auto r = results.begin(); r != results.end(); r++) {
            std::cout << strand << " ";
            r->print();
        }
    }

}

int main(int argc, char** argv) {
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

    std::cerr << "Building forward FMI\n";
    NanoFMI fwd_fmi(model.kmer_count(), fwd_ids, tally_gap);

    std::cerr << "Building reverse FMI\n";
    NanoFMI rev_fmi(model.kmer_count(), rev_ids, tally_gap);

    /* navigate through all the read files provided */
    for (int f = 4; f < argc; f++) {
        std::vector<Event> events;
        if (not fast5::File::is_valid_file(argv[f])) {
            std::cerr << "<" << argv[f] << "> is not a valid file. skipping... " << std::endl;
            continue;
        }
        fast5::File file;
        try
        {
            file.open(argv[f]);
            assert(file.is_open());

            if (file.have_eventdetection_events()) {
                events = file.get_eventdetection_events();
            } else {
                std::cerr << "file " << argv[f] << " does not contain events. skipping..." << std::endl;
                continue;
            }

            NormParams scale = model.get_norm_params(events);

            Timer timer;
            
            std::cout << "== " << argv[f] << " ==\n";

            align_kmers("rev", model, rev_fmi, events, scale);

            std::cerr << "Reverse time: " << timer.lap() << "\n";

            align_kmers("fwd", model, fwd_fmi, events, scale);

            std::cout << "== done ==\n";

            std::cerr << "Forward time: " << timer.lap() << "\n";

        }
        catch (hdf5_tools::Exception& e)
        {
            std::cerr << "hdf5 error: " << e.what() << std::endl;
        }
    }
}
