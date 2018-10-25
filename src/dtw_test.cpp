#include <iostream>
#include <math.h>

#include "kmer_model.hpp"
#include "fast5.hpp"
#include "arg_parse.hpp"
#include "timer.hpp"
#include "dtw.hpp"

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

enum Opt {MODEL     = 'm',
          REF_FASTA = 'x',
          FAST5_NAME = 'i',
          DIAG_COEF = 'd',
          VERT_COEF = 'v',
          HORZ_COEF = 'h',
          NORM_SCALE = 'a',
          NORM_SHIFT = 'b',
          ROW_SUBSTR = 'R',
          COL_SUBSTR = 'C',
          EVENTS = 'E',
          ALN_REV = 'V'};

int main(int argc, char** argv) {

    ArgParse args("DTW tester");
    args.add_string(Opt::MODEL, "model", "/home-4/skovaka1@jhu.edu/code/nanopore_aligner/kmer_models/r9.4_180mv_450bps_6mer/5mers_pre.txt", "Nanopore kmer model");
    args.add_string(Opt::REF_FASTA, "ref_fasta", "", "");
    args.add_string(Opt::FAST5_NAME, "fast5_name", "", "");
    args.add_double(Opt::DIAG_COEF, "diag_coef", 1, "");
    args.add_double(Opt::VERT_COEF, "vert_coef", 1, "");
    args.add_double(Opt::HORZ_COEF, "horz_coef", 1, "");
    args.add_double(Opt::NORM_SCALE, "norm_scale", 0, "");
    args.add_double(Opt::NORM_SHIFT, "norm_shift", 0, "");
    args.add_flag(Opt::ROW_SUBSTR, "row_substr", "");
    args.add_flag(Opt::COL_SUBSTR, "col_substr", "");
    args.add_flag(Opt::EVENTS, "create_events", "");
    args.add_flag(Opt::ALN_REV, "aln_rev", "");
    args.parse_args(argc, argv);

    const KmerModel model(args.get_string(Opt::MODEL), false);

    EventParams params = event_detection_defaults;
    params.threshold1 = 1.4;
    params.threshold2 = 1.1;
    EventDetector event_detector(params);

    SubSeqDTW substr = SubSeqDTW::NONE;
    if (args.get_flag(Opt::ROW_SUBSTR) && args.get_flag(Opt::COL_SUBSTR)) {
        std::cerr << "Error: only use one substr\n";
        return 1;
    }

    if (args.get_flag(Opt::ROW_SUBSTR)) {
        substr = SubSeqDTW::ROW;
    } else if (args.get_flag(Opt::COL_SUBSTR)) {
        substr = SubSeqDTW::COL;
    }

    auto costfn = [&model](float e, u16 k) {return -model.event_match_prob(e,k);};

    std::ifstream ref_in(args.get_string(Opt::REF_FASTA));
    std::vector<u16> fwd_kmers, rev_kmers;
    model.parse_fasta(ref_in, fwd_kmers, rev_kmers);

    //std::cerr << "Opening fast5\n";

    fast5::File fast5_file;
    if (!open_fast5(args.get_string(Opt::FAST5_NAME), fast5_file)) {
        return 1;
    }

    //std::cerr << "Reading raw\n";
    std::vector<float> signal = fast5_file.get_raw_samples();
    if (args.get_flag(Opt::EVENTS)) {
        std::cerr << "Detecting events\n";
        auto events = event_detector.get_all_events(signal);
        signal.clear();
        for (auto e : events) {
            signal.push_back(e.mean);
        }
    }

    NormParams norm = {args.get_double(Opt::NORM_SHIFT),
                       args.get_double(Opt::NORM_SCALE)};
    if (norm.scale <= 0) {
        norm = model.get_norm_params(signal);
    }
    //std::vector<float> signal_unnorm(signal);
    model.normalize(signal, norm);

    //Takes up too much space :(
    if (signal.size() > 500000) {
        return 1;
    }

    //std::cerr << "Aligning fwd\n";

    DTW<float, u16, decltype(costfn)> fwd_dtw(signal, fwd_kmers, 
                                            costfn, substr, 
                                            args.get_double(Opt::DIAG_COEF),
                                            args.get_double(Opt::VERT_COEF),
                                            args.get_double(Opt::HORZ_COEF));
    if (args.get_flag(Opt::ALN_REV)) {
        //std::cerr << "Aligning rev\n";
        DTW<float, u16, decltype(costfn)> rev_dtw(signal, rev_kmers, 
                                                costfn, substr, 
                                                args.get_double(Opt::DIAG_COEF),
                                                args.get_double(Opt::VERT_COEF),
                                                args.get_double(Opt::HORZ_COEF));
        std::cout << fwd_dtw.mean_score() << "\t" << rev_dtw.mean_score() << "\n";
        if (fwd_dtw.mean_score() > rev_dtw.mean_score()) {
            //rev_dtw.print_path();
        }
            return 0;
    }

    std::cout << fwd_dtw.mean_score() << "\t"
              << signal.size() << "\t" 
              << fwd_kmers.size() << "\t"
              << norm.scale << "\t"
              << norm.shift << "\n";
    fwd_dtw.print_path();
    return 0;
}
