/* sigalign.cpp 
 * align oxford nanopore event-level signals directly to the reference without
 * basecalling */

#include <iostream>
#include <string>
#include <unordered_map>
#include <math.h>
#include <unistd.h>

#include "fast5.hpp"
#include "seed_graph.hpp"
#include "seed_tracker.hpp"
#include "kmer_model.hpp"
#include "nano_fmi.hpp"
#include "timer.h"

const char *OPTIONS = "m:r:t:l:e:d:y:s:f:";

const char MODEL_OPT      = 'm',
           REFERENCE_OPT  = 'r',
           TALLY_DIST_OPT = 't',
           READ_LIST_OPT  = 'l',
           EVENT_PT_OPT   = 'e',
           SEED_PT_OPT    = 'd',
           STAY_PT_OPT    = 'y',
           SEED_LEN_OPT   = 's',
           STAY_FRAC_OPT  = 'f';

//TODO: relative paths
const std::string MODEL_DEF = "/home-4/skovaka1@jhu.edu/code/nanopore_aligner/kmer_models/r9.2_180mv_250bps_6mer/template_median68pA.model";

const int TALLY_DIST_DEF = 1,
          SEED_LEN_DEF   = 32;

const double EVENT_PT_DEF  = -5.29,
             SEED_PT_DEF   = -3.35,
             STAY_PT_DEF   = -8.343,
             STAY_FRAC_DEF = 0.7;

bool get_events(std::string filename, std::vector<Event> &events);

int main(int argc, char** argv) {
    std::string model_filename = MODEL_DEF,
                ref_filename = "",
                reads_filename = "";

    int tally_dist = TALLY_DIST_DEF,
        seed_len   = SEED_LEN_DEF;  

    double event_prob_thr = EVENT_PT_DEF,
           seed_prob_thr  = SEED_PT_DEF,
           stay_prob_thr  = STAY_PT_DEF,
           stay_frac      = STAY_FRAC_DEF;

    char o;
    while ( (o = getopt(argc, argv, OPTIONS)) != -1 ) {
        switch (o) {
        case MODEL_OPT:
            model_filename = std::string(optarg);
            break;
        case REFERENCE_OPT:
            ref_filename = std::string(optarg);
            break;
        case TALLY_DIST_OPT:
            tally_dist = atoi(optarg);
            break;
        case READ_LIST_OPT:   
            reads_filename = std::string(optarg);
            break;
        case EVENT_PT_OPT:
            event_prob_thr = atof(optarg);
            break;
        case SEED_PT_OPT: 
            seed_prob_thr = atof(optarg);
            break;
        case STAY_PT_OPT:  
            stay_prob_thr = atof(optarg);
            break;
        case SEED_LEN_OPT:  
            seed_len = atoi(optarg);
            break;
        case STAY_FRAC_OPT: 
            stay_frac = atof(optarg);
            break;
        default:
            std::cerr << "Error: unrecognized option '" << o << "'\n";
            return 1;
        }
    }

    std::cerr << "Loading model\n";
    KmerModel model(model_filename);

    std::cerr << "Parsing fasta\n";
    std::vector<mer_id> fwd_ids, rev_ids;
    std::ifstream ref_file(ref_filename);
    model.parse_fasta(ref_file, fwd_ids, rev_ids);

    std::cerr << "Building forward FMI\n";
    NanoFMI fwd_fmi(model.kmer_count(), fwd_ids, tally_dist);

    std::cerr << "Building reverse FMI\n";
    NanoFMI rev_fmi(model.kmer_count(), rev_ids, tally_dist);

    AlnParams aln_params = {seed_len,
                            event_prob_thr,
                            seed_prob_thr,
                            stay_prob_thr,
                            stay_frac};

    SeedGraph rev_sg(model, rev_fmi, aln_params),
              fwd_sg(model, fwd_fmi, aln_params);

    std::ifstream reads_file(reads_filename);

    if (!reads_file) {
        std::cerr << "Error: couldn't open '" << reads_filename << "'\n";
        return 1;
    }

    std::string read_filename;

    while (getline(reads_file, read_filename)) {
        std::vector<Event> events;

        
        if (!get_events(read_filename, events)) {
            continue;
        }

        NormParams norm = model.get_norm_params(events);

        fwd_sg.new_read(events.size(), norm);
        rev_sg.new_read(events.size(), norm);

        SeedTracker fwd_tracker, rev_tracker;

        Timer timer;
        
        std::cout << "== " << read_filename << " ==\n";
        std::cout << "== " << events.size() << " events ==\n";

        std::cerr << read_filename << "\n";

        int status_step = events.size() / 10,
            status = 0;

        for (int i = events.size()-1; i >= seed_len; i--) {
            std::vector<Result> fwd_seeds = fwd_sg.add_event(events[i]),
                                rev_seeds = rev_sg.add_event(events[i]);
            fwd_tracker.add_seeds(fwd_seeds);
            rev_tracker.add_seeds(rev_seeds);

            if (status == status_step) {
                int prog = (int) ((100.0 * (events.size() - i)) / events.size());
                std::cerr << prog << "%  (" << timer.get() / 1000 << ")\n";
                status = 0;
            } else {
                status++;
            }
        }
        std::cerr << "100% (" << timer.get() / 1000 << ")\n";


        std::vector<ReadAln> fwd_alns = fwd_tracker.get_alignments(1),
                             rev_alns = rev_tracker.get_alignments(1);

        for (int i = 0; i < fwd_alns.size(); i++) {
            std::cout << "fwd ";
            fwd_alns[i].print();
        }

        for (int i = 0; i < rev_alns.size(); i++) {
            std::cout << "rev ";
            rev_alns[i].print();
        }

        std::cout << "== " << timer.lap() / 1000 << " sec ==\n";
    }
}

bool get_events(std::string filename, std::vector<Event> &events) {
    events.clear();

    if (!fast5::File::is_valid_file(filename)) {
        std::cerr << "Error: '" << filename << "' is not a valid file \n";
    }

    try {
        fast5::File file;
        file.open(filename);
        
        if (!file.is_open()) {  
            std::cerr << "Error: unable to open '" 
                      << filename << "'\n";

        } else if (!file.have_eventdetection_events()) {
            std::cerr << "Error: file '" << filename
                      << "' does not contain events\n";

        } else {
            events = file.get_eventdetection_events();
            return true;
        }

    } catch (hdf5_tools::Exception& e) {
        std::cerr << "Error: hdf5 exception '" << e.what() << "'\n";
    }

    return false;
}
