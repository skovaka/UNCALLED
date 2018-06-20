/* sigalign.cpp 
 * align oxford nanopore event-level signals directly to the reference without
 * basecalling */

#include <iostream>
#include <string>
#include <unordered_map>
#include <math.h>
#include <unistd.h>
#include <algorithm>
#include <sdsl/suffix_arrays.hpp>
#include <iomanip>

#include "fast5.hpp"
#include "arg_parse.hpp"
#include "graph_aligner.hpp"
#include "forest_aligner.hpp"
#include "leaf_aligner.hpp"
#include "leaf_arr_aligner.hpp"
#include "seed_tracker.hpp"  
#include "range.hpp"
#include "kmer_model.hpp"
#include "sdsl_fmi.hpp"
#include "base_fmi.hpp"
#include "timer.hpp"

extern "C" {
    #include "event_detection/event_detection.h"
    #include "event_detection/fast5_interface.h"
}

#define GRAPH_ALN 0
#define FOREST_ALN 1
#define LEAF_ALN 2
#define LEAF_ARR_ALN 3

#ifndef ALN_TYPE
#define ALN_TYPE LEAF_ALN
#endif

#define BASE_FMI 0
#define SDSL_FMI 1

#ifndef FMI_TYPE
#define FMI_TYPE SDSL_FMI
#endif

std::string FWD_STR = "fwd",
            REV_STR = "rev";

//TODO: relative paths
const std::string MODEL_DEF = "/home-4/skovaka1@jhu.edu/code/nanopore_aligner/kmer_models/r9.2_180mv_250bps_6mer/template_median68pA.model";

bool get_events(std::ostream &err, std::string filename, std::vector<Event> &events);

enum Opt {
    MODEL       = 'm',
    REFERENCE   = 'f',
    INDEX_PREFIX   = 'x',
    READ_LIST   = 'i',
    OUT_PREFIX     = 'o',

    TALLY_DIST  = 't',

    MIN_ALN_LEN     = 'a',
    MIN_ALN_CONF    = 'u', //u for until
    PATH_WIN_LEN    = 'w',
    MIN_REP_LEN     = 'r',
    MAX_REP_COPY    = 'c',
    MAX_PATHS       = 'p',

    STAY_FRAC       = 's',
    MAX_CONSEC_STAY = 'y',
    MAX_IGNORES     = 'g',
    MAX_SKIPS       = 'k',

    EVENT_PROBS = 'E',
    WINDOW_PROB = 'W'
};

int main(int argc, char** argv) {
    ArgParse args("UNCALLED: Utility for Nanopore Current Alignment to Large Expanses of DNA");

    args.add_string(Opt::MODEL, "model", "/home-4/skovaka1@jhu.edu/code/nanopore_aligner/kmer_models/r9.2_180mv_250bps_6mer/template_median68pA.model", "Nanopore kmer model");
    args.add_string(Opt::REFERENCE,    "reference",    "",     "");
    args.add_string(Opt::READ_LIST,    "read_list",    "",     "");
    args.add_string(Opt::INDEX_PREFIX, "index_prefix",   "",     "");
    args.add_string(Opt::OUT_PREFIX,   "out_prefix",   "./",     "");

    args.add_int(   Opt::TALLY_DIST, "tally_dist",   128,     "");

    args.add_int(   Opt::MIN_ALN_LEN,  "min_aln_len", 25, "");
    args.add_double(Opt::MIN_ALN_CONF, "min_aln_conf", 2, "");
    args.add_int(   Opt::PATH_WIN_LEN, "path_window_len", 22, "");
    args.add_int(   Opt::MIN_REP_LEN,  "min_repeat_len", 0, "");
    args.add_int(   Opt::MAX_REP_COPY, "max_repeat_copy", 100, "");
    args.add_int(   Opt::MAX_PATHS,    "max_paths", 10000, "");

    args.add_double(Opt::STAY_FRAC,       "stay_frac",    0.5,    "");
    args.add_int(   Opt::MAX_CONSEC_STAY, "max_consec_stay", 8, "");
    args.add_int(   Opt::MAX_IGNORES,     "max_ignores",  0,    "");
    args.add_int(   Opt::MAX_SKIPS,       "max_skips",    0,    "");

    args.add_string(Opt::EVENT_PROBS, "event_probs", "-2.25_100-2.4_5-4.0_1-10.0", "");
    args.add_double(Opt::WINDOW_PROB, "window_prob", -3.75, "");

    args.parse_args(argc, argv);

    std::string prefix = args.get_string(Opt::OUT_PREFIX) + args.get_string(Opt::EVENT_PROBS) + args.get_param_str();
    std::ofstream seeds_out(prefix + "_seeds.txt");
    std::ofstream alns_out(prefix + "_aln.txt");
    std::ofstream time_out(prefix + "_time.txt");

    std::cerr << "Writing:" << "\n"
              << prefix << "_seeds.txt" << "\n"
              << prefix << "_aln.txt" << "\n"
              << prefix << "_time.txt" << "\n";

    KmerModel model(args.get_string(Opt::MODEL));

    AlnParams aln_params(model,
                         args.get_int(Opt::PATH_WIN_LEN),
                         args.get_int(Opt::MIN_REP_LEN),
                         args.get_int(Opt::MAX_REP_COPY),
                         args.get_int(Opt::MAX_PATHS),
                         (float) args.get_double(Opt::STAY_FRAC),
                         args.get_int(Opt::MAX_CONSEC_STAY),
                         args.get_int(Opt::MAX_IGNORES),
                         args.get_int(Opt::MAX_SKIPS),
                         args.get_string(Opt::EVENT_PROBS),
                         (float) args.get_double(Opt::WINDOW_PROB));
    
    #if FMI_TYPE == BASE_FMI
    BaseFMI fwd_fmi, rev_fmi;
    #else
    SdslFMI fwd_fmi, rev_fmi;
    #endif

    if (!args.get_string(Opt::INDEX_PREFIX).empty()) {
        std::string p = args.get_string(Opt::INDEX_PREFIX);

        #if FMI_TYPE == BASE_FMI

        std::ifstream fwd_in(p + "fwdFM.txt"),
                      rev_in(p + "revFM.txt");

        std::cerr << "Reading forward FMI\n";
        fwd_fmi = BaseFMI(fwd_in, args.get_int(Opt::TALLY_DIST));

        std::cerr << "Reading reverse FMI\n";
        rev_fmi = BaseFMI(rev_in, args.get_int(Opt::TALLY_DIST));

        #else

        std::cerr << "Reading forward FMI\n";
        fwd_fmi = SdslFMI(p + "fwdFM.idx");

        std::cerr << "Reading reverse FMI\n";
        rev_fmi = SdslFMI(p + "revFM.idx");

        #endif

    } else {
        std::cerr << "Parsing fasta\n";
        std::ifstream ref_file(args.get_string(Opt::REFERENCE));
        std::string fwd_bases, rev_bases;
        parse_fasta(ref_file, fwd_bases, rev_bases, false);

        #if FMI_TYPE == BASE_FMI
        std::cerr << "Building forward FMI\n";
        fwd_fmi = BaseFMI(fwd_bases, args.get_int(Opt::TALLY_DIST));

        std::cerr << "Building reverse FMI\n";
        rev_fmi = BaseFMI(rev_bases, args.get_int(Opt::TALLY_DIST));

        #else

        std::cerr << "Building forward FMI\n";
        fwd_fmi.construct(fwd_bases);

        std::cerr << "Building reverse FMI\n";
        rev_fmi.construct(rev_bases);

        #endif
    }


    #if ALN_TYPE == GRAPH_ALN
    GraphAligner
    #elif ALN_TYPE == FOREST_ALN
    ForestAligner
    #elif ALN_TYPE == LEAF_ALN
    LeafAligner
    #elif ALN_TYPE == LEAF_ARR_ALN
    LeafArrAligner
    #endif
        rev_sg(rev_fmi, aln_params, "rev"),
        fwd_sg(fwd_fmi, aln_params, "fwd");

    std::ifstream reads_file(args.get_string(Opt::READ_LIST));

    if (!reads_file) {
        std::cerr << "Error: couldn't open '" 
                  << args.get_string(Opt::READ_LIST) << "'\n";
        return 1;
    }

    unsigned int min_aln_len = args.get_int(Opt::MIN_ALN_LEN);
    float min_aln_conf = (float) args.get_double(Opt::MIN_ALN_CONF);
    bool read_until = min_aln_conf > 0;

    std::string read_line, read_filename;

    std::cerr << read_until << " Aligning...\n";

    while (getline(reads_file, read_line)) {

        unsigned int aln_st = 0, aln_en = 0;

        unsigned int i = read_line.find("\t");

        if (i == std::string::npos) {
            read_filename = read_line;
        } else {
            read_filename = read_line.substr(0, i);

            unsigned int j = read_line.find("\t", i+1);

            if (j == std::string::npos) {
                aln_st = atoi(read_line.substr(i+1).c_str());
            } else {
                aln_st = atoi(read_line.substr(i+1, j).c_str());
                aln_en = atoi(read_line.substr(j+1).c_str());
            }
        }
            

        std::vector<Event> events;
        
        if (!get_events(std::cerr, read_filename, events)) {
            continue;
        }

        NormParams norm = model.get_norm_params(events);
        model.normalize_events(events, norm);

        if (aln_en == aln_st) {
            aln_en = events.size() - 1;
        }

        //NormParams norm = model.get_norm_params(events);


        unsigned int aln_len = aln_en - aln_st + 1;

        fwd_sg.new_read(aln_len);
        rev_sg.new_read(aln_len);

        SeedTracker fwd_tracker(aln_len), rev_tracker(aln_len);

        Timer read_timer;

        #ifdef VERBOSE_TIME
        //Timer event_timer;

        time_out << std::fixed << std::setprecision(3);
        #else
        unsigned int status_step = aln_len / 10,
                     status = 0;
        #endif
        
        seeds_out << "== " << read_filename << " ==\n";
        alns_out << "== " << read_filename << " ==\n";
        alns_out << "== scale/shift: " << norm.scale << " " << norm.shift << " ==\n";
        alns_out << "== aligning " 
                 << aln_st << "-" 
                 << aln_en << " of " 
                 << events.size() << " events ==\n";

        #ifndef VERBOSE_TIME
        time_out << "== " << read_filename << " ==\n";
        #endif

        float **all_probs = new float*[aln_len];

        unsigned int e = aln_en;

        for (; e >= aln_st && e <= aln_en; e--) {

            all_probs[e-aln_st] = new float[model.kmer_count()];
            for (Kmer kmer = 0; kmer < model.kmer_count(); kmer++) {
                all_probs[e-aln_st][kmer] = model.event_match_prob(events[e], kmer);
            }

            //#ifdef VERBOSE_TIME
            //time_out << e;
            ////event_timer.reset();
            //#endif

            std::vector<Result> fwd_seeds = fwd_sg.add_event(all_probs[e-aln_st], 
                                                             seeds_out,
                                                             time_out);


            std::vector<Result> rev_seeds = rev_sg.add_event(all_probs[e-aln_st], 
                                                            seeds_out,
                                                            time_out);

            //#ifdef VERBOSE_TIME
            //event_timer.reset();
            //#endif

            ReadAln fa = fwd_tracker.add_seeds(fwd_seeds);

            //#ifdef VERBOSE_TIME
            //time_out << std::setw(23) << event_timer.lap();
            //#endif

            ReadAln ra = rev_tracker.add_seeds(rev_seeds);

            if (aln_en - e > 50000) {
                break;
            }

            //#ifdef VERBOSE_TIME
            //time_out << std::setw(23) << event_timer.lap();
            //#endif

            if (read_until && std::max(fa.total_len_, ra.total_len_) > min_aln_len) {
                ReadAln a = fa;
                if (fa.total_len_ < ra.total_len_) {
                    a = ra;
                }

                if(fwd_tracker.check_ratio(a, min_aln_conf) &&
                   rev_tracker.check_ratio(a, min_aln_conf)) {

                    //#ifdef VERBOSE_TIME
                    //time_out << std::setw(23) << event_timer.lap() << std::endl;
                    //#endif

                    break;
                }
            }

            #ifndef VERBOSE_TIME
            if (status == status_step) {
                unsigned int prog = (unsigned int) ((100.0 * (aln_len - e)) / aln_len);
                time_out << prog << "%  (" << read_timer.get() / 1000 << ")\n";
                status = 0;
                time_out.flush();
                seeds_out.flush();
                alns_out.flush();
            } else {
                status++;
            }
            #endif

            time_out.flush();
        }

        #ifdef VERBOSE_TIME
        //time_out << std::setw(23) << event_timer.lap() << std::endl;
        #if ALN_TYPE == LEAF_ARR_ALN
        time_out << std::setw(10) << (fwd_sg.loop1_time_ + rev_sg.loop1_time_) 
                 << std::setw(10) << (fwd_sg.sort_time_ + rev_sg.sort_time_) 
                 << std::setw(10) << (fwd_sg.loop2_time_ + rev_sg.loop2_time_) 
                 << std::setw(10) << (fwd_sg.fullsource_time_ + rev_sg.fullsource_time_) 
                 << std::setw(10) << (fwd_sg.fmrs_time_ + rev_sg.fmrs_time_) 
                 << std::setw(10) << (fwd_sg.fmsa_time_ + rev_sg.fmsa_time_) << "\n";
        #endif

        #else
        time_out << "100% (" << read_timer.get() / 1000 << ")\n";

        time_out << "== END ==\n";
        #endif

        fwd_tracker.print(alns_out, FWD_STR, (size_t) 5);
        rev_tracker.print(alns_out, REV_STR, (size_t) 5);

        std::cout.flush();
        
        time_out.flush();
        seeds_out.flush();
        alns_out.flush();

        alns_out  << "== " << read_timer.lap() << " ms, " 
                  << (aln_en - e + 1) << " events ==\n";

        seeds_out << "== " << read_timer.lap() / 1000 << " sec ==\n";

        for (e = e-aln_st+1; e < aln_len; e++) {
            delete[] all_probs[e];
        }
        delete[] all_probs;
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
