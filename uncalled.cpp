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
#include "leaf_arr_aligner.hpp"
#include "seed_tracker.hpp"  
#include "range.hpp"
#include "kmer_model.hpp"
#include "bwa_fmi.hpp"
#include "sdsl_fmi.hpp"
#include "base_fmi.hpp"
#include "timer.hpp"

#define GRAPH_ALN 0
#define FOREST_ALN 1
#define LEAF_ALN 2
#define LEAF_ARR_ALN 3

#ifndef ALN_TYPE
#define ALN_TYPE LEAF_ALN
#endif

#define BASE_FMI 0
#define SDSL_FMI 1
#define BWA_FMI 2

#ifndef FMI_TYPE
#define FMI_TYPE BWA_FMI
#endif

std::string FWD_STR = "fwd",
            REV_STR = "rev";

//TODO: relative paths
const std::string MODEL_DEF = "/home-4/skovaka1@jhu.edu/code/nanopore_aligner/kmer_models/r9.2_180mv_250bps_6mer/template_median68pA.model";

bool get_events(std::string filename, std::vector<Event> &events, EventDetector &ed);

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

    KmerModel model(args.get_string(Opt::MODEL), true);

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
    BaseFMI fmi;
    #elif FMI_TYPE == SDSL_FMI
    SdslFMI fmi;
    #else
    BwaFMI fmi;
    #endif

    std::cerr << "Reading FMI\n";

    #if FMI_TYPE == BASE_FMI
    std::ifstream fmi_in(args.get_string(Opt::INDEX_PREFIX));
    fmi = BaseFMI(fmi_in, args.get_int(Opt::TALLY_DIST));

    #elif FMI_TYPE == SDSL_FMI
    fmi = SdslFMI(args.get_string(Opt::INDEX_PREFIX));

    #else
    fmi = BwaFMI(args.get_string(Opt::INDEX_PREFIX));
    #endif

    #if ALN_TYPE == GRAPH_ALN
    GraphAligner
    #elif ALN_TYPE == FOREST_ALN
    ForestAligner
    #elif ALN_TYPE == LEAF_ALN
    LeafAligner
    #elif ALN_TYPE == LEAF_ARR_ALN
    LeafArrAligner
    #endif
        graph(fmi, aln_params);

    std::ifstream reads_file(args.get_string(Opt::READ_LIST));

    if (!reads_file) {
        std::cerr << "Error: couldn't open '" 
                  << args.get_string(Opt::READ_LIST) << "'\n";
        return 1;
    }

    unsigned int min_aln_len = args.get_int(Opt::MIN_ALN_LEN);
    float min_aln_conf = (float) args.get_double(Opt::MIN_ALN_CONF);
    bool read_until = min_aln_conf > 0;

    detector_param params = event_detection_defaults;
    params.threshold1 = 1.4;
    params.threshold2 = 1.1;
    EventDetector event_detector(params, 30, 150);

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
        
        //TODO: get_raw_sample_dataset
        if (!get_events(read_filename, events, event_detector)) {
            continue;
        }

        NormParams norm = model.get_norm_params(events);
        model.normalize_events(events, norm);

        if (aln_en == aln_st) {
            aln_en = events.size() - 1;
        }

        unsigned int aln_len = aln_en - aln_st + 1;

        graph.new_read(aln_len);

        SeedTracker tracker(aln_len);

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


        u32 e = aln_st;
        std::vector<float> kmer_probs(model.kmer_count());

        for (; e >= aln_st && e <= aln_en; e++) {

            for (u16 kmer = 0; kmer < model.kmer_count(); kmer++) {
                kmer_probs[kmer] = model.event_match_prob(events[e], kmer);
            }

            std::vector<Result> seeds = graph.add_event(kmer_probs, 
                                                        seeds_out,
                                                        time_out);


            ReadAln new_aln = tracker.add_seeds(seeds);

            //if (aln_en - e > 20000) {
            //    break;
            //}

            if (read_until && 
                new_aln.total_len_ > min_aln_len && 
                tracker.check_ratio(new_aln, min_aln_conf)) {
                    break;
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
        time_out << std::setw(10) << (graph.loop1_time_) 
                 << std::setw(10) << (graph.sort_time_) 
                 << std::setw(10) << (graph.loop2_time_)
                 << std::setw(10) << (graph.fullsource_time_)
                 << std::setw(10) << (graph.fmrs_time_)
                 << std::setw(10) << (graph.fmsa_time_) << "\n";
        #endif

        #else
        time_out << "100% (" << read_timer.get() / 1000 << ")\n";

        time_out << "== END ==\n";
        #endif

        tracker.print(alns_out, FWD_STR, (size_t) 10);

        std::cout.flush();
        
        time_out.flush();
        seeds_out.flush();
        alns_out.flush();

        alns_out  << "== " << read_timer.get() << " ms, " 
                  << (e - aln_st + 1) << " events ==\n";

        seeds_out << "== " << read_timer.get() << " ms ==\n";
    }

}

bool get_events(std::string filename, std::vector<Event> &events, EventDetector &ed) {
    if (!fast5::File::is_valid_file(filename)) {
        std::cerr << "Error: '" << filename << "' is not a valid file \n";
    }

    try {
        fast5::File file;
        file.open(filename);
        
        if (!file.is_open()) {  
            std::cerr << "Error: unable to open '" 
                      << filename << "'\n";
            return false;
        }

        events = ed.get_all_events(file.get_raw_samples());
        return true;
        
        /*else if (!file.have_eventdetection_events()) {

            std::cerr << "Error: no events\n";

        } else {
            events = file.get_eventdetection_events();
            return true;
        }*/

    } catch (hdf5_tools::Exception& e) {
        std::cerr << "Error: hdf5 exception '" << e.what() << "'\n";
        return false;
    }

    return false;
}
