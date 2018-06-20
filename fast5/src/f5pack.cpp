//
// Part of: https://github.com/mateidavid/fast5
//
// Copyright (c) 2015-2017 Matei David, Ontario Institute for Cancer Research
// MIT License
//

#include <cassert>
#include <iostream>
#include <iomanip>
#include <string>

#include <tclap/CmdLine.h>
#include "logger.hpp"

#include "fast5.hpp"
#include "File_Packer.hpp"

using namespace std;

namespace opts
{
    using namespace TCLAP;
    string description = "Pack an ONT fast5 file.";
    CmdLine cmd_parser(description);
    //
    MultiArg< string > log_level("", "log", "Log level. (default: info)", false, "string", cmd_parser);
    MultiSwitchArg extra_verbosity("v", "", "Increase verbosity", cmd_parser);
    //
    SwitchArg al_drop("", "al-drop", "Drop basecall alignment data.", cmd_parser);
    SwitchArg al_copy("", "al-copy", "Copy basecall alignment data.", cmd_parser);
    SwitchArg al_unpack("", "al-unpack", "Unpack basecall alignment data.", cmd_parser);
    SwitchArg al_pack("", "al-pack", "Pack basecall alignment data.", cmd_parser);
    //
    SwitchArg ev_drop("", "ev-drop", "Drop basecall event data.", cmd_parser);
    SwitchArg ev_copy("", "ev-copy", "Copy basecall event data.", cmd_parser);
    SwitchArg ev_unpack("", "ev-unpack", "Unpack basecall event data.", cmd_parser);
    SwitchArg ev_pack("", "ev-pack", "Pack basecall event data.", cmd_parser);
    //
    SwitchArg fq_drop("", "fq-drop", "Drop basecall fastq data.", cmd_parser);
    SwitchArg fq_copy("", "fq-copy", "Copy basecall fastq data.", cmd_parser);
    SwitchArg fq_unpack("", "fq-unpack", "Unpack basecall fatsq data.", cmd_parser);
    SwitchArg fq_pack("", "fq-pack", "Pack basecall fastq data.", cmd_parser);
    //
    SwitchArg ed_drop("", "ed-drop", "Drop event detection data.", cmd_parser);
    SwitchArg ed_copy("", "ed-copy", "Copy event detection data.", cmd_parser);
    SwitchArg ed_unpack("", "ed-unpack", "Unpack event detection data.", cmd_parser);
    SwitchArg ed_pack("", "ed-pack", "Pack event detection data.", cmd_parser);
    //
    SwitchArg rw_drop("", "rw-drop", "Drop raw samples data.", cmd_parser);
    SwitchArg rw_copy("", "rw-copy", "Copy raw samples data.", cmd_parser);
    SwitchArg rw_unpack("", "rw-unpack", "Unpack raw samples data.", cmd_parser);
    SwitchArg rw_pack("", "rw-pack", "Pack raw samples data.", cmd_parser);
    //
    ValueArg< unsigned > p_model_state_bits("", "p-model-state-bits", "P_Model_State bits to keep.", false, fast5::File_Packer::default_p_model_state_bits(), "int", cmd_parser);
    ValueArg< unsigned > qv_bits("", "qv-bits", "QV bits to keep.", false, fast5::File_Packer::max_qv_bits(), "int", cmd_parser);
    SwitchArg no_check("n", "no-check", "Don't check packing.", cmd_parser);
    SwitchArg force("f", "force", "Overwrite output file if it exists.", cmd_parser);
    //
    SwitchArg fastq("", "fastq", "Pack fastq data, drop rest.", cmd_parser);
    SwitchArg archive("", "archive", "Pack raw saples data, drop rest.", cmd_parser);
    SwitchArg unpack("u", "unpack", "Unpack files.", cmd_parser);
    SwitchArg pack("p", "pack", "Pack files (default, if no other pack/unpack/copy options).", cmd_parser);
    //
    UnlabeledValueArg< string > input_fn("input", "Input fast5 file.", true, "", "file", cmd_parser);
    UnlabeledValueArg< string > output_fn("output", "Output fast5 file.", true, "", "file", cmd_parser);
} // opts


int main(int argc, char * argv[])
{
    opts::cmd_parser.parse(argc, argv);
    // set log levels
    auto default_level = (int)logger::level::info + opts::extra_verbosity.getValue();
    logger::Logger::set_default_level(default_level);
    logger::Logger::set_levels_from_options(opts::log_level, &clog);
    // print options
    LOG(info) << "program: " << opts::cmd_parser.getProgramName() << endl;
    LOG(info) << "version: " << opts::cmd_parser.getVersion() << endl;
    LOG(info) << "args: " << opts::cmd_parser.getOrigArgv() << endl;
    // what to pack/unpack
    if (opts::pack + opts::unpack + opts::archive + opts::fastq > 1)
    {
        LOG_EXIT << "at most one of --pack/--unpack/--archive/--fastq may be given" << endl;
    }
    if (opts::rw_pack + opts::rw_unpack + opts::rw_copy + opts::rw_drop > 1)
    {
        LOG_EXIT << "at most one of --rw-pack/--rw-unpack/--rw-copy/--rw-drop may be given" << endl;
    }
    if (opts::ed_pack + opts::ed_unpack + opts::ed_copy + opts::ed_drop > 1)
    {
        LOG_EXIT << "at most one of --ed-pack/--ed-unpack/--ed-copy/--ed-drop may be given" << endl;
    }
    if (opts::fq_pack + opts::fq_unpack + opts::fq_copy + opts::fq_drop > 1)
    {
        LOG_EXIT << "at most one of --fq-pack/--fq-unpack/--fq-copy/--fq-drop may be given" << endl;
    }
    if (opts::ev_pack + opts::ev_unpack + opts::ev_copy + opts::ev_drop > 1)
    {
        LOG_EXIT << "at most one of --ev-pack/--ev-unpack/--ev-copy/--ev-drop may be given" << endl;
    }
    if (opts::al_pack + opts::al_unpack + opts::al_copy + opts::al_drop > 1)
    {
        LOG_EXIT << "at most one of --al-pack/--al-unpack/--al-copy/--al-drop may be given" << endl;
    }
    if (opts::pack + opts::unpack + opts::archive + opts::fastq
        + opts::rw_pack + opts::rw_unpack + opts::rw_copy + opts::rw_drop
        + opts::ed_pack + opts::ed_unpack + opts::ed_copy + opts::ed_drop
        + opts::fq_pack + opts::fq_unpack + opts::fq_copy + opts::fq_drop
        + opts::ev_pack + opts::ev_unpack + opts::ev_copy + opts::ev_drop
        + opts::al_pack + opts::al_unpack + opts::al_copy + opts::al_drop
        == 0)
    {
        opts::pack.set(true);
    }
    if (opts::pack)
    {
        opts::rw_pack.set(true);
        opts::ed_pack.set(true);
        opts::fq_pack.set(true);
        opts::ev_pack.set(true);
        opts::al_pack.set(true);
    }
    else if (opts::unpack)
    {
        opts::rw_unpack.set(true);
        opts::ed_unpack.set(true);
        opts::fq_unpack.set(true);
        opts::ev_unpack.set(true);
        opts::al_unpack.set(true);
    }
    if (opts::archive)
    {
        opts::rw_pack.set(true);
    }
    if (opts::fastq)
    {
        opts::fq_pack.set(true);
    }
    if (opts::rw_pack + opts::rw_unpack + opts::rw_copy + opts::rw_drop == 0) opts::rw_drop.set(true);
    if (opts::ed_pack + opts::ed_unpack + opts::ed_copy + opts::ed_drop == 0) opts::ed_drop.set(true);
    if (opts::fq_pack + opts::fq_unpack + opts::fq_copy + opts::fq_drop == 0) opts::fq_drop.set(true);
    if (opts::ev_pack + opts::ev_unpack + opts::ev_copy + opts::ev_drop == 0) opts::ev_drop.set(true);
    if (opts::al_pack + opts::al_unpack + opts::al_copy + opts::al_drop == 0) opts::al_drop.set(true);
    LOG(info) << "rw: " << (opts::rw_pack? "pack" : opts::rw_unpack? "unpack" : opts::rw_copy? "copy" : "drop") << endl;
    LOG(info) << "ed: " << (opts::ed_pack? "pack" : opts::ed_unpack? "unpack" : opts::ed_copy? "copy" : "drop") << endl;
    LOG(info) << "fq: " << (opts::fq_pack? "pack" : opts::fq_unpack? "unpack" : opts::fq_copy? "copy" : "drop") << endl;
    LOG(info) << "ev: " << (opts::ev_pack? "pack" : opts::ev_unpack? "unpack" : opts::ev_copy? "copy" : "drop") << endl;
    LOG(info) << "al: " << (opts::al_pack? "pack" : opts::al_unpack? "unpack" : opts::al_copy? "copy" : "drop") << endl;
    LOG(info) << "check: " << (not opts::no_check? "yes" : "no") << endl;
    // set File_Packer options
    int rw_policy = (opts::rw_pack? 1 : opts::rw_unpack? 2 : opts::rw_copy? 3 : 0);
    int ed_policy = (opts::ed_pack? 1 : opts::ed_unpack? 2 : opts::ed_copy? 3 : 0);
    int fq_policy = (opts::fq_pack? 1 : opts::fq_unpack? 2 : opts::fq_copy? 3 : 0);
    int ev_policy = (opts::ev_pack? 1 : opts::ev_unpack? 2 : opts::ev_copy? 3 : 0);
    int al_policy = (opts::al_pack? 1 : opts::al_unpack? 2 : opts::al_copy? 3 : 0);
    fast5::File_Packer fp(rw_policy, ed_policy, fq_policy, ev_policy, al_policy);
    fp.set_check(not opts::no_check);
    fp.set_force(opts::force);
    fp.set_qv_bits(opts::qv_bits);
    fp.set_p_model_state_bits(opts::p_model_state_bits);
    fp.run(opts::input_fn, opts::output_fn);
    auto cnt = fp.get_counts();
    cout
        << std::fixed << std::setprecision(2)
        << "bp_seq_count\t" << cnt.bp_seq_count << "\n"
        << "rs_count\t" << cnt.rs_count << "\n"
        << "rs_bits\t" << (double)cnt.rs_bits/cnt.rs_count << "\n"
        << "ed_count\t" << cnt.ed_count << "\n"
        << "ed_skip_bits\t" << (double)cnt.ed_skip_bits/cnt.ed_count << "\n"
        << "ed_len_bits\t" << (double)cnt.ed_len_bits/cnt.ed_count << "\n"
        << "fq_count\t" << cnt.fq_count << "\n"
        << "fq_bp_bits\t" << (double)cnt.fq_bp_bits/cnt.fq_count << "\n"
        << "fq_qv_bits\t" << (double)cnt.fq_qv_bits/cnt.fq_count << "\n"
        << "ev_count\t" << cnt.ev_count << "\n"
        << "ev_rel_skip_bits\t" << (double)cnt.ev_rel_skip_bits/cnt.ev_count << "\n"
        << "ev_skip_bits\t" << (double)cnt.ev_skip_bits/cnt.ev_count << "\n"
        << "ev_len_bits\t" << (double)cnt.ev_len_bits/cnt.ev_count << "\n"
        << "ev_move_bits\t" << (double)cnt.ev_move_bits/cnt.ev_count << "\n"
        << "ev_p_model_state_bits\t" << (double)cnt.ev_p_model_state_bits/cnt.ev_count << "\n"
        << "al_count\t" << cnt.al_count << "\n"
        << "al_template_step_bits\t" << (double)cnt.al_template_step_bits/cnt.al_count << "\n"
        << "al_complement_step_bits\t" << (double)cnt.al_complement_step_bits/cnt.al_count << "\n"
        << "al_move_bits\t" << (double)cnt.al_move_bits/cnt.al_count << "\n";
}
