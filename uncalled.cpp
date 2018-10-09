/* sigalign.cpp 
 * align oxford nanopore event-level signals directly to the reference without
 * basecalling */

#include <iostream>
#include <string>
#include <unordered_map>
#include <math.h>
#include <unistd.h>
#include <algorithm>
#include <iomanip>

#include "fast5.hpp"
#include "arg_parse.hpp"
#include "leaf_arr_aligner.hpp"
#include "seed_tracker.hpp"  
#include "range.hpp"
#include "kmer_model.hpp"
#include "bwa_fmi.hpp"
#include "base_fmi.hpp"
#include "timer.hpp"

std::string FWD_STR = "fwd",
            REV_STR = "rev";

//TODO: relative paths
const std::string MODEL_DEF = "/home-4/skovaka1@jhu.edu/code/nanopore_aligner/kmer_models/r9.2_180mv_250bps_6mer/template_median68pA.model";

bool open_fast5(const std::string &filename, fast5::File &file);

enum Opt {
    MODEL       = 'm',
    REFERENCE   = 'f',
    INDEX_PREFIX   = 'x',
    READ_LIST   = 'i',
    OUT_PREFIX     = 'o',

    TALLY_DIST  = 't',

    MIN_ALN_LEN     = 'a',
    PATH_WIN_LEN    = 'w',
    MIN_REP_LEN     = 'r',
    MAX_REP_COPY    = 'c',
    MAX_PATHS       = 'p',

    MIN_MEAN_CONF   = 'M', //u for until
    MIN_TOP_CONF    = 'T', //u for until

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
    args.add_int(   Opt::PATH_WIN_LEN, "path_window_len", 22, "");
    args.add_int(   Opt::MIN_REP_LEN,  "min_repeat_len", 0, "");
    args.add_int(   Opt::MAX_REP_COPY, "max_repeat_copy", 100, "");
    args.add_int(   Opt::MAX_PATHS,    "max_paths", 10000, "");

    args.add_double(Opt::MIN_MEAN_CONF, "min_mean_conf", 6.67, "");
    args.add_double(Opt::MIN_TOP_CONF,  "min_top_conf", 2, "");

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
    

    std::cerr << "Reading FMI\n";
    BwaFMI fmi(args.get_string(Opt::INDEX_PREFIX));

    LeafArrAligner graph(fmi, aln_params);

    std::ifstream reads_file(args.get_string(Opt::READ_LIST));

    if (!reads_file) {
        std::cerr << "Error: couldn't open '" 
                  << args.get_string(Opt::READ_LIST) << "'\n";
        return 1;
    }

    unsigned int min_aln_len = args.get_int(Opt::MIN_ALN_LEN);
    float min_mean_conf = (float) args.get_double(Opt::MIN_MEAN_CONF),
          min_top_conf = (float) args.get_double(Opt::MIN_TOP_CONF);
    bool read_until = min_mean_conf > 0 || min_top_conf > 0;
    std::cout << min_mean_conf << " "
              << min_top_conf << " "
              << read_until << std::endl;

    detector_param params = event_detection_defaults;
    params.threshold1 = 1.4;
    params.threshold2 = 1.1;
    EventDetector event_detector(params, 30, 150);

    std::string read_line, read_filename;

    std::cerr << read_until << " Aligning...\n";

    std::cout << fmi.bns_->l_pac << " len\n";

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
            
        fast5::File fast5_file;
        if (!open_fast5(read_filename, fast5_file)) {
            continue;
        }

        //TODO: get_raw_sample_dataset
        auto fast5_info = fast5_file.get_raw_samples_params();
        auto raw_samples = fast5_file.get_raw_samples();
        auto events = event_detector.get_all_events(raw_samples);

        NormParams norm = model.get_norm_params(events);
        model.normalize(events, norm);

        if (aln_en == aln_st) {
            aln_en = events.size() - 1;
        }

        unsigned int aln_len = aln_en - aln_st + 1;

        graph.new_read(aln_len);
        

        SeedTracker tracker(fmi.size(), 
                            min_mean_conf, 
                            min_top_conf, 
                            min_aln_len, 
                            args.get_int(Opt::PATH_WIN_LEN));


        Timer read_timer;

        #ifdef VERBOSE_TIME
        Timer tracker_timer;
        double tracker_time = 0;
        time_out << std::fixed << std::setprecision(3);
        #else
        unsigned int status_step = aln_len / 10,
                     status = 0;
        time_out << "== " << read_filename << " ==\n";
        #endif
        
        seeds_out << "== " << read_filename << " ==\n";


        u32 e = aln_st;
        std::vector<float> kmer_probs(model.kmer_count());

        ReadAln final_aln;

        for (; e >= aln_st && e <= aln_en; e++) {

            for (u16 kmer = 0; kmer < model.kmer_count(); kmer++) {
                kmer_probs[kmer] = model.event_match_prob(events[e], kmer);
            }

            std::vector<Result> seeds = graph.add_event(kmer_probs, 
                                                        seeds_out,
                                                        time_out);

            #ifdef VERBOSE_TIME
            tracker_timer.reset();
            #endif

            final_aln = tracker.add_seeds(seeds);

            #ifdef VERBOSE_TIME
            tracker_time += tracker_timer.get();
            #endif

            //if (aln_en - e > 20000) {
            //    break;
            //}

            if (read_until && final_aln.is_valid()) {
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

        double aln_time = read_timer.get();

        #ifdef VERBOSE_TIME
            //time_out << std::setw(23) << event_timer.lap() << std::endl;
            #if ALN_TYPE == LEAF_ARR_ALN
            time_out << std::setw(10) << (graph.loop1_time_) 
                     << std::setw(10) << (graph.sort_time_) 
                     << std::setw(10) << (graph.loop2_time_)
                     << std::setw(10) << (graph.fullsource_time_)
                     << std::setw(10) << (graph.fmrs_time_)
                     << std::setw(10) << (graph.fmsa_time_)
                     << std::setw(10) << tracker_time << "\n";
            #endif
        #else
        time_out << "100% (" << aln_time / 1000 << ")\n";

        time_out << "== END ==\n";
        #endif

        if (read_until) {
            if (final_aln.is_valid()) {

                u64 st, en;
                std::string strand;
                if (final_aln.ref_st_ < fmi.size() / 2) {
                    st = final_aln.ref_st_;
                    en = final_aln.ref_en_.end_;
                    strand = "rev";
                } else {
                    st = fmi.size() - final_aln.ref_en_.end_ - 1;
                    en = fmi.size() - final_aln.ref_st_ - 1;
                    strand = "fwd";
                }

                i32 rid = bns_pos2rid(fmi.bns_, st);
                std::string ref_name(fmi.bns_->anns[rid].name);
                st -= fmi.bns_->anns[rid].offset;
                en -= fmi.bns_->anns[rid].offset;

                alns_out << "aligned\t" 
                         << final_aln.evt_st_ << "-" << final_aln.evt_en_  << "\t"
                         << ref_name << "\t" << strand << "\t" 
                         << st << "-" << en << "\t"
                         << final_aln.total_len_ << "\t";
            } else {
                alns_out << "unaligned\t0-" << aln_len << "\tNULL\tNULL\t0-0\t0\t";
            }
            alns_out 
                     << aln_time << "\t" 
                     << aln_len << "\t"
                     << fast5_info.read_id << "\t"
                     << read_filename << std::endl;
        } else {
            tracker.print(alns_out, 10);
            alns_out << aln_len << " events processed" << std::endl;
        }

        seeds_out << "== " << aln_time << " ms ==\n";

        time_out.flush();
        seeds_out.flush();
        alns_out.flush();
    }

}

bool open_fast5(const std::string &filename, fast5::File &file) {
    if (!fast5::File::is_valid_file(filename)) {
        std::cerr << "Error: '" << filename << "' is not a valid file \n";
    }

    try {
        file.open(filename);
        
        if (!file.is_open()) {  
            std::cerr << "Error: unable to open '" << filename << "'\n";
            return false;
        }

        return true;
        
    } catch (hdf5_tools::Exception& e) {
        std::cerr << "Error: hdf5 exception '" << e.what() << "'\n";
        return false;
    }

    return false;
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
