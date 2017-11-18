/* sigalign.cpp 
 * align oxford nanopore event-level signals directly to the reference without
 * basecalling */

#include <iostream>
#include <string>
#include <unordered_map>
#include <math.h>
#include <unistd.h>

#include "fast5.hpp"
#include "arg_parse.hpp"
#include "seed_graph.hpp"
#include "seed_tracker.hpp"
#include "kmer_model.hpp"
#include "kmer_fmi.hpp"
#include "timer.h"


std::string FWD_STR = "fwd",
            REV_STR = "rev";

EventAligner::EventAligner(std::string model_fname,
                           std::stirng ref_fname,
                           int min_seed_len,
                           int anchor_len,
                           int max_ignores,
                           int max_skips,
                           int stay_frac,
                           int anchor_evpr,
                           int extend_evpr,
                           int seed_pr,
                           int stay_pr) {
                           

    err_out << "Loading model\n";
    KmerModel model(args.get_string(Opt::MODEL));

    err_out << "Parsing fasta\n";
    std::vector<mer_id> fwd_ids, rev_ids;
    std::ifstream ref_file(args.get_string(Opt::REFERENCE));
    model.parse_fasta(ref_file, fwd_ids, rev_ids);

    err_out << "Building forward FMI\n";
    KmerFMI fwd_fmi(model.kmer_count(), fwd_ids, args.get_int(Opt::TALLY_DIST));

    err_out << "Building reverse FMI\n";
    KmerFMI rev_fmi(model.kmer_count(), rev_ids, args.get_int(Opt::TALLY_DIST));

    AlnParams aln_params(model,
                         args.get_int(Opt::MIN_SEED_LEN),
                         args.get_int(Opt::ANCHOR_LEN),
                         args.get_int(Opt::MAX_IGNORES),
                         args.get_int(Opt::MAX_SKIPS),
                         args.get_double(Opt::STAY_FRAC),
                         args.get_double(Opt::ANCHOR_EVPR),
                         args.get_double(Opt::EXTEND_EVPR),
                         args.get_double(Opt::SEED_PR),
                         args.get_double(Opt::STAY_PR));

    SeedGraph rev_sg(rev_fmi, aln_params, "rev"),
              fwd_sg(fwd_fmi, aln_params, "fwd");


    std::ifstream reads_file(args.get_string(Opt::READ_LIST));

    if (!reads_file) {
        err_out << "Error: couldn't open '" 
                  << args.get_string(Opt::READ_LIST) << "'\n";
        return 1;
    }

    std::string read_line, read_filename;


    while (getline(reads_file, read_line)) {

        int aln_st = 0, aln_en = -1;

        int i = read_line.find("\t");

        if (i == std::string::npos) {
            read_filename = read_line;
        } else {
            read_filename = read_line.substr(0, i);

            int j = read_line.find("\t", i+1);

            if (j == std::string::npos) {
                aln_st = atoi(read_line.substr(i+1).c_str());
            } else {
                aln_st = atoi(read_line.substr(i+1, j).c_str());
                aln_en = atoi(read_line.substr(j+1).c_str());
            }
        }
            

        std::vector<Event> events;
        
        if (!get_events(err_out, read_filename, events)) {
            continue;
        }

        if (aln_en < aln_st) {
            aln_en = events.size() - 1;
        }

        NormParams norm = model.get_norm_params(events);

        err_out << "NORM: " << norm.scale << " " << norm.shift << "\n";

        int aln_len = aln_en - aln_st + 1;

        fwd_sg.new_read(aln_len, norm);
        rev_sg.new_read(aln_len, norm);

        SeedTracker fwd_tracker, rev_tracker;

        Timer timer;
        
        seeds_out << "== " << read_filename << " ==\n";
        alns_out << "== " << read_filename << " ==\n";
        alns_out << "== aligning " 
                 << aln_st << "-" 
                 << aln_en << " of " 
                 << events.size() << " events ==\n";

        err_out << read_filename << "\n";

        int status_step = aln_len / 10,
            status = 0;

        for (int i = aln_en; i >= aln_st; i--) {
            std::vector<Result> fwd_seeds = fwd_sg.add_event(events[i], seeds_out),
                                rev_seeds = rev_sg.add_event(events[i], seeds_out);
            fwd_tracker.add_seeds(fwd_seeds);
            rev_tracker.add_seeds(rev_seeds);

            //alns_out << events[i].mean << " " << events[i].stdv << "\n";

            if (status == status_step) {
                int prog = (int) ((100.0 * (aln_len - i)) / aln_len);
                err_out << prog << "%  (" << timer.get() / 1000 << ")\n";
                status = 0;
                err_out.flush();
                seeds_out.flush();
                alns_out.flush();
            } else {
                status++;
            }
        }
        err_out << "100% (" << timer.get() / 1000 << ")\n";

        fwd_tracker.print(alns_out, FWD_STR, 5);
        rev_tracker.print(alns_out, REV_STR, 5);

        err_out.flush();
        seeds_out.flush();
        alns_out.flush();

        alns_out << "== " << timer.lap() / 1000 << " sec ==\n";
        seeds_out << "== " << timer.lap() / 1000 << " sec ==\n";
    }
}

bool get_events(std::ostream &err, std::string filename, std::vector<Event> &events) {
    events.clear();

    if (!fast5::File::is_valid_file(filename)) {
        err << "Error: '" << filename << "' is not a valid file \n";
    }

    try {
        fast5::File file;
        file.open(filename);
        
        if (!file.is_open()) {  
            err << "Error: unable to open '" 
                      << filename << "'\n";

        } else if (!file.have_eventdetection_events()) {
            err << "Error: file '" << filename
                      << "' does not contain events\n";

        } else {
            events = file.get_eventdetection_events();
            return true;
        }

    } catch (hdf5_tools::Exception& e) {
        err << "Error: hdf5 exception '" << e.what() << "'\n";
    }

    return false;
}
