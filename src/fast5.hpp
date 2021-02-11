//
// Part of: https://github.com/mateidavid/fast5
//
// Copyright (c) 2015-2017 Matei David, Ontario Institute for Cancer Research
// MIT License
//

#ifndef __FAST5_HPP
#define __FAST5_HPP

#include <algorithm>
#include <cassert>
#include <cmath>
#include <exception>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <array>
#include <set>
#include <map>
#include <stdexcept>

#include "fast5/logger.hpp"
#include "fast5/fast5_version.hpp"
#include "fast5/hdf5_tools.hpp"
#include "fast5/Huffman_Packer.hpp"
#include "fast5/Bit_Packer.hpp"

#define MAX_K_LEN 8

namespace
{
    inline static std::string array_to_string(std::array< char, MAX_K_LEN > const & a)
    {
        return std::string(a.begin(), std::find(a.begin(), a.end(), '\0'));
    }
}

namespace fast5
{

typedef hdf5_tools::File::Attr_Map Attr_Map;

struct Channel_Id_Params
{
    std::string channel_number;
    double digitisation;
    double offset;
    double range;
    double sampling_rate;
    Channel_Id_Params()
        : channel_number(""),
          digitisation(0.0),
          offset(0.0),
          range(0.0),
          sampling_rate(0.0) {}
    void read(hdf5_tools::File const & f, std::string const & p)
    {
        f.read(p + "/channel_number", channel_number);
        f.read(p + "/digitisation", digitisation);
        f.read(p + "/offset", offset);
        f.read(p + "/range", range);
        f.read(p + "/sampling_rate", sampling_rate);
    }
    void write(hdf5_tools::File const & f, std::string const & p) const
    {
        f.write_attribute(p + "/channel_number", channel_number);
        f.write_attribute(p + "/digitisation", digitisation);
        f.write_attribute(p + "/offset", offset);
        f.write_attribute(p + "/range", range);
        f.write_attribute(p + "/sampling_rate", sampling_rate);
    }
}; // struct Channel_Id_Params

typedef Attr_Map Tracking_Id_Params;

typedef Attr_Map Sequences_Params;

typedef Attr_Map Context_Tags_Params;

typedef float Raw_Sample;
typedef int16_t Raw_Int_Sample;

struct Raw_Samples_Params
{
    std::string read_id;
    long long read_number;
    long long start_mux;
    long long start_time;
    long long duration;
    friend bool operator == (Raw_Samples_Params const & lhs, Raw_Samples_Params const & rhs)
    {
        return (lhs.read_id == rhs.read_id
                and lhs.read_number == rhs.read_number
                and lhs.start_mux == rhs.start_mux
                and lhs.start_time == rhs.start_time
                and lhs.duration == rhs.duration);
    }
    void read(hdf5_tools::File const & f, std::string const & p)
    {
        f.read(p + "/read_id", read_id);
        f.read(p + "/read_number", read_number);
        f.read(p + "/start_mux", start_mux);
        f.read(p + "/start_time", start_time);
        f.read(p + "/duration", duration);
    }
    void write(hdf5_tools::File const & f, std::string const & p) const
    {
        f.write_attribute(p + "/read_id", read_id);
        f.write_attribute(p + "/read_number", read_number);
        f.write_attribute(p + "/start_mux", start_mux);
        f.write_attribute(p + "/start_time", start_time);
        f.write_attribute(p + "/duration", duration);
    }
}; // struct Raw_Samples_Params

typedef std::pair< std::vector< Raw_Int_Sample >, Raw_Samples_Params > Raw_Int_Samples_Dataset;
typedef std::pair< std::vector< Raw_Sample >, Raw_Samples_Params > Raw_Samples_Dataset;

struct Raw_Samples_Pack
{
    Huffman_Packer::Code_Type signal;
    Attr_Map signal_params;
    //
    Raw_Samples_Params params;
    //
    void read(hdf5_tools::File const & f, std::string const & p)
    {
        f.read(p + "/Signal", signal);
        signal_params = f.get_attr_map(p + "/Signal");
        params.read(f, p + "/params");
    }
    void write(hdf5_tools::File const & f, std::string const & p) const
    {
        f.write_dataset(p + "/Signal", signal);
        f.add_attr_map(p + "/Signal", signal_params);
        params.write(f, p + "/params");
    }
}; // struct Raw_Samples_Pack

struct EventDetection_Event
{
    double mean;
    double stdv;
    long long start;
    long long length;
    friend bool operator == (EventDetection_Event const & lhs, EventDetection_Event const & rhs)
    {
        return lhs.mean == rhs.mean
            and lhs.stdv == rhs.stdv
            and lhs.start == rhs.start
            and lhs.length == rhs.length;
    }
    static hdf5_tools::Compound_Map const & compound_map()
    {
        static hdf5_tools::Compound_Map m;
        static bool inited = false;
        if (not inited)
        {
            m.add_member("mean", &EventDetection_Event::mean);
            m.add_member("start", &EventDetection_Event::start);
            m.add_member("length", &EventDetection_Event::length);
            m.add_member("stdv", &EventDetection_Event::stdv);
            inited = true;
        }
        return m;
    }
    static hdf5_tools::Compound_Map const & alt_compound_map()
    {
        static hdf5_tools::Compound_Map m;
        static bool inited = false;
        if (not inited)
        {
            m.add_member("mean", &EventDetection_Event::mean);
            m.add_member("start", &EventDetection_Event::start);
            m.add_member("length", &EventDetection_Event::length);
            m.add_member("variance", &EventDetection_Event::stdv);
            inited = true;
        }
        return m;
    }
}; // struct EventDetection_Event

struct EventDetection_Events_Params
{
    std::string read_id;
    long long read_number;
    long long scaling_used;
    long long start_mux;
    long long start_time;
    long long duration;
    double median_before;
    unsigned abasic_found;
    friend bool operator == (EventDetection_Events_Params const & lhs, EventDetection_Events_Params const & rhs)
    {
        return (lhs.read_id == rhs.read_id
                and lhs.read_number == rhs.read_number
                and lhs.scaling_used == rhs.scaling_used
                and lhs.start_mux == rhs.start_mux
                and lhs.start_time == rhs.start_time
                and lhs.duration == rhs.duration
                and ((std::isnan(lhs.median_before) and std::isnan(rhs.median_before))
                     or lhs.median_before == rhs.median_before)
                and lhs.abasic_found == rhs.abasic_found);
    }
    void read(hdf5_tools::File const & f, std::string const & p)
    {
        auto a_v = f.get_attr_list(p);
        std::set< std::string > a_s(a_v.begin(), a_v.end());
        f.read(p + "/read_number", read_number);
        f.read(p + "/scaling_used", scaling_used);
        f.read(p + "/start_mux", start_mux);
        f.read(p + "/start_time", start_time);
        f.read(p + "/duration", duration);
        // optional fields
        if (a_s.count("read_id"))
        {
            f.read(p + "/read_id", read_id);
        }
        if (a_s.count("median_before"))
        {
            f.read(p + "/median_before", median_before);
        }
        else
        {
            median_before = std::nan("");
        }
        if (a_s.count("abasic_found"))
        {
            f.read(p + "/abasic_found", abasic_found);
        }
        else
        {
            abasic_found = 2;
        }
    }
    void write(hdf5_tools::File const & f, std::string const & p) const
    {
        f.write_attribute(p + "/read_number", read_number);
        f.write_attribute(p + "/scaling_used", scaling_used);
        f.write_attribute(p + "/start_mux", start_mux);
        f.write_attribute(p + "/start_time", start_time);
        f.write_attribute(p + "/duration", duration);
        if (not read_id.empty()) f.write_attribute(p + "/read_id", read_id);
        if (not std::isnan(median_before)) f.write_attribute(p + "/median_before", median_before);
        if (abasic_found < 2) f.write_attribute(p + "/abasic_found", abasic_found);
    }
}; // struct EventDetection_Events_Params

typedef std::pair< std::vector< EventDetection_Event >, EventDetection_Events_Params > EventDetection_Events_Dataset;

struct EventDetection_Events_Pack
{
    Huffman_Packer::Code_Type skip;
    Attr_Map skip_params;
    Huffman_Packer::Code_Type len;
    Attr_Map len_params;
    //
    EventDetection_Events_Params params;
    //
    void read(hdf5_tools::File const & f, std::string const & p)
    {
        f.read(p + "/Skip", skip);
        skip_params = f.get_attr_map(p + "/Skip");
        f.read(p + "/Len", len);
        len_params = f.get_attr_map(p + "/Len");
        params.read(f, p + "/params");
    }
    void write(hdf5_tools::File const & f, std::string const & p) const
    {
        f.write_dataset(p + "/Skip", skip);
        f.add_attr_map(p + "/Skip", skip_params);
        f.write_dataset(p + "/Len", len);
        f.add_attr_map(p + "/Len", len_params);
        params.write(f, p + "/params");
    }
}; // struct EventDetection_Events_Pack

//
// This struct represents the expected signal measured
// given the kmer sequence that is in the pore when the
// the observations are made. A pore model consists
// of 1024 of these entries (one per 5-mer) and global
// shift/scaling params.
//
struct Basecall_Model_State
{
    double level_mean;
    double level_stdv;
    double sd_mean;
    double sd_stdv;
    std::array< char, MAX_K_LEN > kmer;
    std::string get_kmer() const { return array_to_string(kmer); }
    friend bool operator == (Basecall_Model_State const & lhs, Basecall_Model_State const & rhs)
    {
        return (lhs.level_mean == rhs.level_mean
                and lhs.level_stdv == rhs.level_stdv
                and lhs.sd_mean == rhs.sd_mean
                and lhs.sd_stdv == rhs.sd_stdv
                and lhs.kmer == rhs.kmer);
    }
    static hdf5_tools::Compound_Map const & compound_map()
    {
        static hdf5_tools::Compound_Map m;
        static bool inited = false;
        if (not inited)
        {
            m.add_member("level_mean", &Basecall_Model_State::level_mean);
            m.add_member("level_stdv", &Basecall_Model_State::level_stdv);
            m.add_member("sd_mean", &Basecall_Model_State::sd_mean);
            m.add_member("sd_stdv", &Basecall_Model_State::sd_stdv);
            m.add_member("kmer", &Basecall_Model_State::kmer);
            inited = true;
        }
        return m;
    }
}; // struct Basecall_Model_State

//
// This struct represents the global transformations
// that must be applied to each Basecall_Model_State
//
struct Basecall_Model_Params
{
    double scale;
    double shift;
    double drift;
    double var;
    double scale_sd;
    double var_sd;
    void read(hdf5_tools::File const & f, std::string const & p)
    {
        f.read(p + "/scale", scale);
        f.read(p + "/shift", shift);
        f.read(p + "/drift", drift);
        f.read(p + "/var", var);
        f.read(p + "/scale_sd", scale_sd);
        f.read(p + "/var_sd", var_sd);
    }
    void write(hdf5_tools::File const & f, std::string const & p) const
    {
        f.write_attribute(p + "/scale", scale);
        f.write_attribute(p + "/shift", shift);
        f.write_attribute(p + "/drift", drift);
        f.write_attribute(p + "/var", var);
        f.write_attribute(p + "/scale_sd", scale_sd);
        f.write_attribute(p + "/var_sd", var_sd);
    }
}; // struct Basecall_Model_Params

struct Basecall_Fastq_Pack
{
    Huffman_Packer::Code_Type bp;
    Attr_Map bp_params;
    Huffman_Packer::Code_Type qv;
    Attr_Map qv_params;
    std::string read_name;
    std::uint8_t qv_bits;
    //
    void read(hdf5_tools::File const & f, std::string const & p)
    {
        f.read(p + "/BP", bp);
        bp_params = f.get_attr_map(p + "/BP");
        f.read(p + "/QV", qv);
        qv_params = f.get_attr_map(p + "/QV");
        f.read(p + "/read_name", read_name);
        f.read(p + "/qv_bits", qv_bits);
    }
    void write(hdf5_tools::File const & f, std::string const & p) const
    {
        f.write_dataset(p + "/BP", bp);
        f.add_attr_map(p + "/BP", bp_params);
        f.write_dataset(p + "/QV", qv);
        f.add_attr_map(p + "/QV", qv_params);
        f.write_attribute(p + "/read_name", read_name);
        f.write_attribute(p + "/qv_bits", qv_bits);
    }
}; // struct Basecall_Fastq_Pack

//
// This struct represents an observed event.
// The members of the struct are the same as
// the fields encoded in the FAST5 file.
//
struct Basecall_Event
{
    double mean;
    double stdv;
    double start;
    double length;
    double p_model_state;
    long long move;
    std::array< char, MAX_K_LEN > model_state;
    std::string get_model_state() const { return array_to_string(model_state); }
    friend bool operator == (Basecall_Event const & lhs, Basecall_Event const & rhs)
    {
        return (lhs.mean == rhs.mean
                and lhs.stdv == rhs.stdv
                and lhs.start == rhs.start
                and lhs.length == rhs.length
                and lhs.p_model_state == rhs.p_model_state
                and lhs.move == rhs.move
                and lhs.model_state == rhs.model_state);
    }
    static hdf5_tools::Compound_Map const & compound_map()
    {
        static hdf5_tools::Compound_Map m;
        static bool inited = false;
        if (not inited)
        {
            m.add_member("mean", &Basecall_Event::mean);
            m.add_member("stdv", &Basecall_Event::stdv);
            m.add_member("start", &Basecall_Event::start);
            m.add_member("length", &Basecall_Event::length);
            m.add_member("p_model_state", &Basecall_Event::p_model_state);
            m.add_member("move", &Basecall_Event::move);
            m.add_member("model_state", &Basecall_Event::model_state);
            inited = true;
        }
        return m;
    }
}; // struct Basecall_Event

struct Basecall_Events_Params
{
    double start_time;
    double duration;
    friend bool operator == (Basecall_Events_Params const & lhs, Basecall_Events_Params const & rhs)
    {
        return (lhs.start_time == rhs.start_time
                and lhs.duration == rhs.duration);
    }
    void read(hdf5_tools::File const & f, std::string const & p)
    {
        if (f.attribute_exists(p + "/start_time"))
        {
            f.read(p + "/start_time", start_time);
        }
        else
        {
            start_time = 0.0;
        }
        if (f.attribute_exists(p + "/duration"))
        {
            f.read(p + "/duration", duration);
        }
        else
        {
            duration = 0.0;
        }
    }
    void write(hdf5_tools::File const & f, std::string const & p) const
    {
        if (start_time > 0.0) f.write_attribute(p + "/start_time", start_time);
        if (duration > 0.0) f.write_attribute(p + "/duration", duration);
    }
};

typedef std::pair< std::vector< Basecall_Event >, Basecall_Events_Params > Basecall_Events_Dataset;

struct Basecall_Events_Pack
{
    Huffman_Packer::Code_Type rel_skip;
    Attr_Map rel_skip_params;
    Huffman_Packer::Code_Type skip;
    Attr_Map skip_params;
    Huffman_Packer::Code_Type len;
    Attr_Map len_params;
    Huffman_Packer::Code_Type move;
    Attr_Map move_params;
    Bit_Packer::Code_Type p_model_state;
    Attr_Map p_model_state_params;
    //
    std::string name;
    std::string version;
    std::string ed_gr;
    long long start_time;
    unsigned state_size;
    double median_sd_temp;
    unsigned p_model_state_bits;
    //
    Basecall_Events_Params params;
    //
    void read(hdf5_tools::File const & f, std::string const & p)
    {
        if (f.dataset_exists(p + "/Rel_Skip"))
        {
            f.read(p + "/Rel_Skip", rel_skip);
            rel_skip_params = f.get_attr_map(p + "/Rel_Skip");
        }
        else
        {
            f.read(p + "/Skip", skip);
            skip_params = f.get_attr_map(p + "/Skip");
            f.read(p + "/Len", len);
            len_params = f.get_attr_map(p + "/Len");
        }
        f.read(p + "/Move", move);
        move_params = f.get_attr_map(p + "/Move");
        f.read(p + "/P_Model_State", p_model_state);
        p_model_state_params = f.get_attr_map(p + "/P_Model_State");
        f.read(p + "/name", name);
        f.read(p + "/version", version);
        f.read(p + "/ed_gr", ed_gr);
        f.read(p + "/start_time", start_time);
        f.read(p + "/state_size", state_size);
        f.read(p + "/median_sd_temp", median_sd_temp);
        f.read(p + "/p_model_state_bits", p_model_state_bits);
        params.read(f, p + "/params");
    }
    void write(hdf5_tools::File const & f, std::string const & p) const
    {
        if (not rel_skip.empty())
        {
            f.write_dataset(p + "/Rel_Skip", rel_skip);
            f.add_attr_map(p + "/Rel_Skip", rel_skip_params);
        }
        else
        {
            f.write_dataset(p + "/Skip", skip);
            f.add_attr_map(p + "/Skip", skip_params);
            f.write_dataset(p + "/Len", len);
            f.add_attr_map(p + "/Len", len_params);
        }
        f.write_dataset(p + "/Move", move);
        f.add_attr_map(p + "/Move", move_params);
        f.write_dataset(p + "/P_Model_State", p_model_state);
        f.add_attr_map(p + "/P_Model_State", p_model_state_params);
        f.write_attribute(p + "/name", name);
        f.write_attribute(p + "/version", version);
        f.write_attribute(p + "/ed_gr", ed_gr);
        f.write_attribute(p + "/start_time", start_time);
        f.write_attribute(p + "/state_size", state_size);
        f.write_attribute(p + "/median_sd_temp", median_sd_temp);
        f.write_attribute(p + "/p_model_state_bits", p_model_state_bits);
        params.write(f, p + "/params");
    }
}; // struct Basecall_Events_Pack

//
// This struct represents a template-to-complement
// match that is emitted by ONT's 2D basecaller
//
struct Basecall_Alignment_Entry
{
    long long template_index;
    long long complement_index;
    std::array< char, MAX_K_LEN > kmer;
    std::string get_kmer() const { return array_to_string(kmer); }
    friend bool operator == (Basecall_Alignment_Entry const & lhs, Basecall_Alignment_Entry const & rhs)
    {
        return lhs.template_index == rhs.template_index
            and lhs.complement_index == rhs.complement_index
            and lhs.kmer == rhs.kmer;
    }
    static hdf5_tools::Compound_Map const & compound_map()
    {
        static hdf5_tools::Compound_Map m;
        static bool inited = false;
        if (not inited)
        {
            m.add_member("template", &Basecall_Alignment_Entry::template_index);
            m.add_member("complement", &Basecall_Alignment_Entry::complement_index);
            m.add_member("kmer", &Basecall_Alignment_Entry::kmer);
            inited = true;
        }
        return m;
    }
}; // struct Basecall_Alignment_Entry

struct Basecall_Alignment_Pack
{
    Bit_Packer::Code_Type template_step;
    Bit_Packer::Code_Params_Type template_step_params;
    Bit_Packer::Code_Type complement_step;
    Bit_Packer::Code_Params_Type complement_step_params;
    Huffman_Packer::Code_Type move;
    Huffman_Packer::Code_Params_Type move_params;
    unsigned template_index_start;
    unsigned complement_index_start;
    unsigned kmer_size;
    //
    void read(hdf5_tools::File const & f, std::string const & p)
    {
        f.read(p + "/Template_Step", template_step);
        template_step_params = f.get_attr_map(p + "/Template_Step");
        f.read(p + "/Complement_Step", complement_step);
        complement_step_params = f.get_attr_map(p + "/Complement_Step");
        f.read(p + "/Move", move);
        move_params = f.get_attr_map(p + "/Move");
        f.read(p + "/template_index_start", template_index_start);
        f.read(p + "/complement_index_start", complement_index_start);
        f.read(p + "/kmer_size", kmer_size);
    }
    void write(hdf5_tools::File const & f, std::string const & p) const
    {
        f.write_dataset(p + "/Template_Step", template_step);
        f.add_attr_map(p + "/Template_Step", template_step_params);
        f.write_dataset(p + "/Complement_Step", complement_step);
        f.add_attr_map(p + "/Complement_Step", complement_step_params);
        f.write_dataset(p + "/Move", move);
        f.add_attr_map(p + "/Move", move_params);
        f.write_attribute(p + "/template_index_start", template_index_start);
        f.write_attribute(p + "/complement_index_start", complement_index_start);
        f.write_attribute(p + "/kmer_size", kmer_size);
    }
};

struct Basecall_Group_Description
{
    std::string name;
    std::string version;
    std::string ed_gr;
    std::string bc_1d_gr;
    bool have_subgroup[3];
    bool have_fastq[3];
    bool have_events[3];
    bool have_model[2];
    bool have_alignment;
    Basecall_Group_Description() :
        have_subgroup{false, false, false},
        have_fastq{false, false, false},
        have_events{false, false, false},
        have_model{false, false},
        have_alignment{false}
    {}
}; // struct Basecall_Group_Description

class File
    : private hdf5_tools::File
{
private:
    typedef hdf5_tools::File Base;
public:
    //
    // Constructors
    //
    File() = default;
    File(std::string const & file_name, bool rw = false) { open(file_name, rw); }

    //
    // Base methods
    //
    using Base::is_open;
    using Base::is_rw;
    using Base::file_name;
    using Base::create;
    using Base::close;
    using Base::get_object_count;
    using Base::is_valid_file;

    //
    // Base method wrappers
    //
    void
    open(std::string const & file_name, bool rw = false)
    {
        Base::open(file_name, rw);
        reload();
    }

    //
    // Access /file_version
    //
    std::string
    file_version() const
    {
        std::string res;
        Base::read(file_version_path(), res);
        return res;
    }

    //
    // Access /UniqueGlobalKey/channel_id
    //
    bool
    have_channel_id_params() const
    {
        return _channel_id_params.sampling_rate > 0.0;
    }
    Channel_Id_Params
    get_channel_id_params() const
    {
        return _channel_id_params;
    }
    void
    add_channel_id_params(Channel_Id_Params const & channel_id_params)
    {
        _channel_id_params = channel_id_params;
        _channel_id_params.write(*this, channel_id_path());
    }
    bool
    have_sampling_rate() const { return have_channel_id_params(); }
    double
    get_sampling_rate() const { return _channel_id_params.sampling_rate; }

    //
    // Access /UniqueGlobalKey/tracking_id
    //
    bool
    have_tracking_id_params() const
    {
        return Base::group_exists(tracking_id_path());
    }
    Tracking_Id_Params
    get_tracking_id_params() const
    {
        return get_attr_map(tracking_id_path());
    }
    void
    add_tracking_id_params(Tracking_Id_Params const & tracking_id_params) const
    {
        add_attr_map(tracking_id_path(), tracking_id_params);
    }
    
    //
    // Access /UniqueGlobalKey/context_tags
    //
    bool
    have_context_tags_params() const
    {
        return Base::group_exists(context_tags_path());
    }

    Context_Tags_Params
    get_context_tags_params() const
    {
        return get_attr_map(context_tags_path());
    }

    void
    add_context_tags_params(Context_Tags_Params const & context_tags_params) const
    {
        add_attr_map(context_tags_path(), context_tags_params);
    }

    //
    // Access /Sequences
    //
    bool
    have_sequences_params() const
    {
        return Base::group_exists(sequences_path());
    }
    Sequences_Params
    get_sequences_params() const
    {
        return get_attr_map(sequences_path());
    }
    void
    add_sequences_params(Sequences_Params const & sequences_params) const
    {
        add_attr_map(sequences_path(), sequences_params);
    }

    //
    // Access Raw Samples
    //
    std::vector< std::string > const &
    get_raw_samples_read_name_list() const
    {
        return _raw_samples_read_names;
    }
    bool
    have_raw_samples(std::string const & rn = std::string()) const
    {
        auto && rn_l = get_raw_samples_read_name_list();
        return (rn.empty()
                ? not rn_l.empty()
                : std::find(rn_l.begin(), rn_l.end(), rn) != rn_l.end());
    }
    bool
    have_raw_samples_unpack(std::string const & rn) const
    {
        return Base::dataset_exists(raw_samples_path(rn));
    }
    bool
    have_raw_samples_pack(std::string const & rn) const
    {
        return Base::group_exists(raw_samples_pack_path(rn));
    }
    Raw_Samples_Params
    get_raw_samples_params(std::string const & rn = std::string()) const
    {
        Raw_Samples_Params res;
        auto && _rn = fill_raw_samples_read_name(rn);
        if (have_raw_samples_unpack(_rn))
        {
            res.read(*this, raw_samples_params_path(_rn));
        }
        else
        {
            res.read(*this, raw_samples_params_pack_path(_rn));
        }
        return res;
    }
    void
    add_raw_samples_params(std::string const & rn, Raw_Samples_Params const & params) const
    {
        std::string p = raw_samples_params_path(rn);
        params.write(*this, p);
    }
    std::vector< Raw_Int_Sample >
    get_raw_int_samples(std::string const & rn = std::string()) const
    {
        std::vector< Raw_Int_Sample > res;
        auto && _rn = fill_raw_samples_read_name(rn);
        if (have_raw_samples_unpack(_rn))
        {
            Base::read(raw_samples_path(_rn), res);
        }
        else if (have_raw_samples_pack(_rn))
        {
            auto rs_pack = get_raw_samples_pack(_rn);
            res = unpack_rw(rs_pack).first;
        }
        return res;
    }
    void
    add_raw_samples(std::string const & rn, std::vector< Raw_Int_Sample > const & rsi)
    {
        Base::write_dataset(raw_samples_path(rn), rsi);
        reload();
    }
    std::vector< Raw_Sample >
    get_raw_samples(std::string const & rn = std::string()) const
    {
        // get raw samples
        auto rsi = get_raw_int_samples(rn);
        // decode levels
        std::vector< Raw_Sample > res;
        res.reserve(rsi.size());
        for (auto int_level : rsi)
        {
            res.push_back(raw_sample_to_float(int_level, _channel_id_params));
        }
        return res;
    }

    //
    // Access EventDetection groups
    //
    std::vector< std::string > const &
    get_eventdetection_group_list() const
    {
        return _eventdetection_groups;
    }
    bool
    have_eventdetection_group(std::string const & gr = std::string()) const
    {
        return (gr.empty()
                ? not _eventdetection_groups.empty()
                : _eventdetection_read_names.count(gr));
    }
    std::vector< std::string > const &
    get_eventdetection_read_name_list(std::string const & gr = std::string()) const
    {
        static const std::vector< std::string > _empty;
        auto && _gr = fill_eventdetection_group(gr);
        return (_eventdetection_read_names.count(_gr)
                ? _eventdetection_read_names.at(_gr)
                : _empty);
    }
    Attr_Map
    get_eventdetection_params(std::string const & gr = std::string()) const
    {
        auto && _gr = fill_eventdetection_group(gr);
        return get_attr_map(eventdetection_group_path(_gr));
    }
    void
    add_eventdetection_params(std::string const & gr, Attr_Map const & am) const
    {
        add_attr_map(eventdetection_group_path(gr), am);
    }

    //
    // Access EventDetection events
    //
    bool
    have_eventdetection_events(
        std::string const & gr = std::string(), std::string const & rn = std::string()) const
    {
        auto && _gr = fill_eventdetection_group(gr);
        auto && _rn = fill_eventdetection_read_name(_gr, rn);
        return (_eventdetection_read_names.count(_gr)
                and std::find(
                    _eventdetection_read_names.at(_gr).begin(),
                    _eventdetection_read_names.at(_gr).end(),
                    _rn)
                != _eventdetection_read_names.at(_gr).end());
    }
    bool
    have_eventdetection_events_unpack(std::string const & gr, std::string const & rn) const
    {
        return Base::dataset_exists(eventdetection_events_path(gr, rn));
    }
    bool
    have_eventdetection_events_pack(std::string const & gr, std::string const & rn) const
    {
        return Base::group_exists(eventdetection_events_pack_path(gr, rn));
    }
    EventDetection_Events_Params
    get_eventdetection_events_params(
        std::string const & gr = std::string(), std::string const & rn = std::string()) const
    {
        EventDetection_Events_Params res;
        auto && _gr = fill_eventdetection_group(gr);
        auto && _rn = fill_eventdetection_read_name(_gr, rn);
        if (have_eventdetection_events_unpack(_gr, _rn))
        {
            res.read(*this, eventdetection_events_params_path(_gr, _rn));
        }
        else if (have_eventdetection_events_pack(_gr, _rn))
        {
            res.read(*this, eventdetection_events_params_pack_path(_gr, _rn));
        }
        return res;
    }
    void
    add_eventdetection_events_params(
        std::string const & gr, std::string const & rn,
        EventDetection_Events_Params const & ede_params) const
    {
        auto p = eventdetection_events_params_path(gr, rn);
        ede_params.write(*this, p);
    }
    std::vector< EventDetection_Event >
    get_eventdetection_events(
        std::string const & gr = std::string(), std::string const & rn = std::string()) const
    {
        std::vector< EventDetection_Event > ede;
        auto && _gr = fill_eventdetection_group(gr);
        auto && _rn = fill_eventdetection_read_name(_gr, rn);
        if (have_eventdetection_events_unpack(_gr, _rn))
        {
            auto p = eventdetection_events_path(_gr, _rn);
            // accept either stdv or variance
            auto meml = get_struct_members(p);
            std::set< std::string > mems(meml.begin(), meml.end());
            if (mems.count("stdv"))
            {
                Base::read(p, ede, EventDetection_Event::compound_map());
            }
            else if (mems.count("variance"))
            {
                Base::read(p, ede, EventDetection_Event::alt_compound_map());
                for (auto & e : ede)
                {
                    e.stdv = std::sqrt(e.stdv);
                }
            }
            else
            {
                LOG_THROW
                    << "neither stdv nor variance found for ed_gr=" << gr;
            }
        }
        else if (have_eventdetection_events_pack(_gr, _rn))
        {
            auto ede_pack = get_eventdetection_events_pack(_gr, _rn);
            if (not have_raw_samples(_rn))
            {
                LOG_THROW_(std::logic_error)
                    << "missing raw samples required to unpack eventdetection events: gr=" << _gr
                    << " rn=" << _rn;
            }
            auto rs_ds = get_raw_samples_dataset(_rn);
            ede = unpack_ed(ede_pack, rs_ds).first;
        }
        return ede;
    } // get_eventdetection_events()
    void
    add_eventdetection_events(
        std::string const & gr, std::string const & rn,
        std::vector< EventDetection_Event > const & ede)
    {
        Base::write_dataset(eventdetection_events_path(gr, rn), ede, EventDetection_Event::compound_map());
        reload();
    }

    //
    // Access Basecall groups
    //
    std::vector< std::string > const &
    get_basecall_group_list() const
    {
        return _basecall_groups;
    }
    bool
    have_basecall_group(std::string const & gr = std::string()) const
    {
        auto && gr_l = get_basecall_group_list();
        return (gr.empty()
                ? not gr_l.empty()
                : std::find(gr_l.begin(), gr_l.end(), gr) != gr_l.end());
    }
    std::vector< std::string > const &
    get_basecall_strand_group_list(unsigned st) const
    {
        return _basecall_strand_groups.at(st);
    }
    bool
    have_basecall_strand_group(unsigned st, std::string const & gr = std::string()) const
    {
        auto && gr_l = get_basecall_strand_group_list(st);
        if (gr.empty())
        {
            return not gr_l.empty();
        }
        if (not _basecall_group_descriptions.count(gr))
        {
            return false;
        }
        else
        {
            return _basecall_group_descriptions.at(gr).have_subgroup[st];
        }
    }
    Basecall_Group_Description const &
    get_basecall_group_description(std::string const & gr) const
    {
        return _basecall_group_descriptions.at(gr);
    }
    std::string const &
    get_basecall_1d_group(std::string const & gr) const
    {
        static std::string const empty;
        return (_basecall_group_descriptions.count(gr)
                ? _basecall_group_descriptions.at(gr).bc_1d_gr
                : empty);
    }
    std::string const &
    get_basecall_eventdetection_group(std::string const & gr) const
    {
        static std::string const empty;
        return (_basecall_group_descriptions.count(gr)
                ? _basecall_group_descriptions.at(gr).ed_gr
                : empty);
    }

    //
    // Access Basecall group params
    //
    Attr_Map
    get_basecall_params(std::string const & gr) const
    {
        return get_attr_map(basecall_group_path(gr));
    }
    void
    add_basecall_params(std::string const & gr, Attr_Map const & am) const
    {
        add_attr_map(basecall_group_path(gr), am);
    }
    //
    // Access Basecall group log
    //
    bool
    have_basecall_log(std::string const & gr) const
    {
        return Base::exists(basecall_log_path(gr));
    }
    std::string
    get_basecall_log(std::string const & gr) const
    {
        std::string res;
        Base::read(basecall_log_path(gr), res);
        return res;
    }
    Attr_Map
    get_basecall_config(std::string const & gr) const
    {
        Attr_Map res;
        if (Base::group_exists(basecall_config_path(gr)))
        {
            res = get_attr_map(basecall_config_path(gr), true);
        }
        return res;
    }
    Attr_Map
    get_basecall_summary(std::string const & gr) const
    {
        Attr_Map res;
        if (Base::group_exists(basecall_summary_path(gr)))
        {
            res = get_attr_map(basecall_summary_path(gr), true);
        }
        return res;
    }

    //
    // Access Basecall fastq
    //
    bool
    have_basecall_fastq(unsigned st, std::string const & gr = std::string()) const
    {
        auto && _gr = fill_basecall_group(st, gr);
        return (_basecall_group_descriptions.count(_gr)
                and _basecall_group_descriptions.at(_gr).have_fastq[st]);
    }
    bool
    have_basecall_fastq_unpack(unsigned st, std::string const & gr) const
    {
        return Base::dataset_exists(basecall_fastq_path(gr, st));
    }
    bool
    have_basecall_fastq_pack(unsigned st, std::string const & gr) const
    {
        return Base::group_exists(basecall_fastq_pack_path(gr, st));
    }
    std::string
    get_basecall_fastq(unsigned st, std::string const & gr = std::string()) const
    {
        std::string res;
        auto && _gr = fill_basecall_group(st, gr);
        if (have_basecall_fastq_unpack(st, _gr))
        {
            Base::read(basecall_fastq_path(_gr, st), res);
        }
        else if (have_basecall_fastq_pack(st, _gr))
        {
            auto fq_pack = get_basecall_fastq_pack(st, _gr);
            res = unpack_fq(fq_pack);
        }
        return res;
    }
    void
    add_basecall_fastq(unsigned st, std::string const & gr, std::string const & fq)
    {
        Base::write(basecall_fastq_path(gr, st), true, fq);
        reload();
    }
    bool
    have_basecall_seq(unsigned st, std::string const & _gr = std::string()) const
    {
        return have_basecall_fastq(st, _gr);
    }
    std::string
    get_basecall_seq(unsigned st, std::string const & _gr = std::string()) const
    {
        return fq2seq(get_basecall_fastq(st, _gr));
    }
    void
    add_basecall_seq(unsigned st, std::string const & gr,
                     std::string const & name, std::string const & seq, int default_qual = 33)
    {
        std::ostringstream oss;
        oss << "@" << name << "\n"
            << seq << "\n"
            << "+\n"
            << std::string(seq.size(), (char)default_qual);
        add_basecall_fastq(st, gr, oss.str());
        reload();
    }

    //
    // Access Basecall model
    //
    bool
    have_basecall_model(unsigned st, std::string const & gr = std::string()) const
    {
        auto && gr_1d = fill_basecall_1d_group(st, gr);
        return (_basecall_group_descriptions.count(gr_1d)
                and _basecall_group_descriptions.at(gr_1d).have_model[st]);
    }
    std::string
    get_basecall_model_file(unsigned st, std::string const & gr = std::string()) const
    {
        std::string res;
        auto && gr_1d = fill_basecall_1d_group(st, gr);
        Base::read(basecall_model_file_path(gr_1d, st), res);
        return res;
    }
    void
    add_basecall_model_file(unsigned st, std::string const & gr, std::string const & file_name) const
    {
        Base::write_attribute(basecall_model_file_path(gr, st), file_name);
    }
    Basecall_Model_Params
    get_basecall_model_params(unsigned st, std::string const & gr = std::string()) const
    {
        Basecall_Model_Params params;
        auto && gr_1d = fill_basecall_1d_group(st, gr);
        std::string path = basecall_model_path(gr_1d, st);
        params.read(*this, path);
        return params;
    }
    void
    add_basecall_model_params(unsigned st, std::string const & gr, Basecall_Model_Params const & params) const
    {
        std::string path = basecall_model_path(gr, st);
        params.write(*this, path);
    }
    std::vector< Basecall_Model_State >
    get_basecall_model(unsigned st, std::string const & gr = std::string()) const
    {
        std::vector< Basecall_Model_State > mod;
        auto && gr_1d = fill_basecall_1d_group(st, gr);
        Base::read(basecall_model_path(gr_1d, st), mod, Basecall_Model_State::compound_map());
        return mod;
    }
    template < typename T >
    void add_basecall_model(unsigned st, std::string const & gr, std::vector< T > const & mod)
    {
        auto && gr_1d = get_basecall_1d_group(gr);
        Base::write_dataset(basecall_model_path(gr_1d, st), mod, Basecall_Model_State::compound_map());
        reload();
    }

    //
    // Access Basecall events
    //
    bool
    have_basecall_events(unsigned st, std::string const & gr = std::string()) const
    {
        auto && gr_1d = fill_basecall_1d_group(st, gr);
        return (_basecall_group_descriptions.count(gr_1d)
                and _basecall_group_descriptions.at(gr_1d).have_events[st]);
    }
    bool
    have_basecall_events_unpack(unsigned st, std::string const & gr) const
    {
        return Base::dataset_exists(basecall_events_path(gr, st));
    }
    bool
    have_basecall_events_pack(unsigned st, std::string const & gr) const
    {
        return Base::group_exists(basecall_events_pack_path(gr, st));
    }
    Basecall_Events_Params
    get_basecall_events_params(unsigned st, std::string const & gr = std::string()) const
    {
        Basecall_Events_Params bce_params;
        auto && gr_1d = fill_basecall_1d_group(st, gr);
        if (have_basecall_events_unpack(st, gr_1d))
        {
            bce_params.read(*this, basecall_events_path(gr_1d, st));
        }
        else if (have_basecall_events_pack(st, gr_1d))
        {
            bce_params.read(*this, basecall_events_params_pack_path(gr_1d, st));
        }
        return bce_params;
    }
    void
    add_basecall_events_params(unsigned st, std::string const & gr,
                               Basecall_Events_Params const & bce_params) const
    {
        auto path = basecall_events_path(gr, st);
        if (not Base::dataset_exists(path))
        {
            LOG_THROW
                << "basecall events must be added before their params";
        }
        bce_params.write(*this, path);
    }
    std::vector< Basecall_Event >
    get_basecall_events(unsigned st, std::string const & gr = std::string()) const
    {
        std::vector< Basecall_Event > res;
        auto && gr_1d = fill_basecall_1d_group(st, gr);
        if (have_basecall_events_unpack(st, gr_1d))
        {
            Base::read(basecall_events_path(gr_1d, st), res, Basecall_Event::compound_map());
        }
        else if (have_basecall_events_pack(st, gr_1d))
        {
            auto ev_pack = get_basecall_events_pack(st, gr_1d);
            if (not have_basecall_seq(st, gr_1d))
            {
                LOG_THROW_(std::logic_error)
                    << "missing fastq required to unpack basecall events: st=" << st
                    << " gr=" << gr_1d;
            }
            auto sq = get_basecall_seq(st, gr_1d);
            if (not ev_pack.ed_gr.empty())
            {
                if (not have_eventdetection_events(ev_pack.ed_gr))
                {
                    LOG_THROW_(std::logic_error)
                        << "missing eventdetection events required to unpack basecall events: st=" << st
                        << " gr=" << gr_1d
                        << " ed_gr=" << ev_pack.ed_gr;
                }
                auto ed = get_eventdetection_events(ev_pack.ed_gr);
                res = unpack_ev(ev_pack, sq, ed, _channel_id_params).first;
            }
            else // ed_gr == "": packed relative to raw samples
            {
                if (not have_raw_samples())
                {
                    LOG_THROW_(std::logic_error)
                        << "missing raw samples required to unpack basecall events: st=" << st
                        << " gr=" << gr_1d;
                }
                auto rs_ds = get_raw_samples_dataset();
                auto ed = unpack_implicit_ed(ev_pack, rs_ds);
                res = unpack_ev(ev_pack, sq, ed, _channel_id_params).first;
            }
        }
        return res;
    }
    template < typename T >
    void
    add_basecall_events(unsigned st, std::string const & gr, std::vector< T > const & ev)
    {
        Base::write_dataset(basecall_events_path(gr, st), ev, T::compound_map());
        reload();
    }

    //
    // Access Basecall alignment
    //
    bool
    have_basecall_alignment(std::string const & gr = std::string()) const
    {
        auto && _gr = fill_basecall_group(2, gr);
        return (_basecall_group_descriptions.count(_gr)
                and _basecall_group_descriptions.at(_gr).have_alignment);
    }
    bool
    have_basecall_alignment_unpack(std::string const & gr) const
    {
        return Base::dataset_exists(basecall_alignment_path(gr));
    }
    bool
    have_basecall_alignment_pack(std::string const & gr) const
    {
        return Base::group_exists(basecall_alignment_pack_path(gr));
    }
    std::vector< Basecall_Alignment_Entry >
    get_basecall_alignment(std::string const & gr = std::string()) const
    {
        std::vector< Basecall_Alignment_Entry > al;
        auto && _gr = fill_basecall_group(2, gr);
        if (have_basecall_alignment_unpack(_gr))
        {
            Base::read(basecall_alignment_path(_gr), al, Basecall_Alignment_Entry::compound_map());
        }
        else if (have_basecall_alignment_pack(_gr)
                 and have_basecall_seq(2, _gr))
        {
            auto al_pack = get_basecall_alignment_pack(_gr);
            auto seq = get_basecall_seq(2, _gr);
            al = unpack_al(al_pack, seq);
        }
        return al;
    }
    void
    add_basecall_alignment(std::string const & gr, std::vector< Basecall_Alignment_Entry > const & al)
    {
        Base::write_dataset(basecall_alignment_path(gr), al, Basecall_Alignment_Entry::compound_map());
        reload();
    }

    //
    // Static helpers
    //
    static inline long long
    time_to_int(double tf, Channel_Id_Params const & cid_params)
    {
        return tf * cid_params.sampling_rate;
    }
    static inline double
    time_to_float(long long ti, Channel_Id_Params const & cid_params)
    {
        return ((long double)ti + .5) / cid_params.sampling_rate;
    }
    static inline float
    raw_sample_to_float(int si, Channel_Id_Params const & cid_params)
    {
        return ((float)si + cid_params.offset)
            * cid_params.range / cid_params.digitisation;
    }
    static std::string
    fq2seq(std::string const & fq)
    {
        return split_fq(fq)[1];
    }
    static std::array< std::string, 4 >
    split_fq(std::string const & fq)
    {
        std::array< std::string, 4 > res = {{"", "", "", ""}};
        size_t i = 0;
        for (int k = 0; k < 4; ++k)
        {
            if (k % 2 == 0) ++i;
            size_t j = fq.find_first_of('\n', i);
            if (j == std::string::npos)
            {
                if (k == 3)
                {
                    j = fq.size();
                }
                else
                {
                    return {{"", "", "", ""}};
                }
            }
            res[k] = fq.substr(i, j - i);
            i = j + 1;
        }
        return res;
    }

private:
    friend struct File_Packer;

    //
    // Cached file data
    //
    Channel_Id_Params _channel_id_params;
    std::vector< std::string > _raw_samples_read_names;
    std::vector< std::string > _eventdetection_groups;
    std::map< std::string, std::vector< std::string > > _eventdetection_read_names;
    std::vector< std::string > _basecall_groups;
    std::map< std::string, Basecall_Group_Description > _basecall_group_descriptions;
    std::array< std::vector< std::string >, 3 > _basecall_strand_groups;

    //
    // Cache updaters
    //
    void
    reload()
    {
        load_channel_id_params();
        load_raw_samples_read_names();
        load_eventdetection_groups();
        load_basecall_groups();
    }
    void
    load_channel_id_params()
    {
        if (not Base::group_exists(channel_id_path())) return;
        _channel_id_params.read(*this, channel_id_path());
    }
    void
    load_raw_samples_read_names()
    {
        _raw_samples_read_names.clear();
        if (not Base::group_exists(raw_samples_root_path())) return;
        auto rn_l = Base::list_group(raw_samples_root_path());
        for (auto const & rn : rn_l)
        {
            if (have_raw_samples_unpack(rn)
                or have_raw_samples_pack(rn))
            {
                _raw_samples_read_names.push_back(rn);
            }
        }
    }
    void
    load_eventdetection_groups()
    {
        _eventdetection_groups.clear();
        _eventdetection_read_names.clear();
        if (not Base::group_exists(eventdetection_root_path())) return;
        auto ed_gr_prefix = eventdetection_group_prefix();
        auto gr_l = Base::list_group(eventdetection_root_path());
        for (auto const & g : gr_l)
        {
            if (g.substr(0, ed_gr_prefix.size()) != ed_gr_prefix) continue;
            std::string gr = g.substr(ed_gr_prefix.size());
            _eventdetection_groups.push_back(gr);
            _eventdetection_read_names[gr] = detect_eventdetection_read_names(gr);
        }
    }
    std::vector< std::string >
    detect_eventdetection_read_names(std::string const & gr) const
    {
        std::vector< std::string > res;
        std::string p = eventdetection_root_path() + "/" + eventdetection_group_prefix() + gr + "/Reads";
        if (not Base::group_exists(p)) return res;
        auto rn_l = Base::list_group(p);
        for (auto const & rn : rn_l)
        {
            if (have_eventdetection_events_unpack(gr, rn)
                or have_eventdetection_events_pack(gr, rn))
            {
                res.push_back(rn);
            }
        }
        return res;
    }
    void
    load_basecall_groups()
    {
        _basecall_groups.clear();
        _basecall_group_descriptions.clear();
        std::for_each(
            _basecall_strand_groups.begin(), _basecall_strand_groups.end(),
            [] (decltype(_basecall_strand_groups)::value_type & v) {
                v.clear();
            });
        if (not Base::group_exists(basecall_root_path())) return;
        auto bc_gr_prefix = basecall_group_prefix();
        auto gr_l = Base::list_group(basecall_root_path());
        for (auto const & g : gr_l)
        {
            if (g.substr(0, bc_gr_prefix.size()) != bc_gr_prefix) continue;
            // found basecall group
            std::string gr = g.substr(bc_gr_prefix.size());
            _basecall_groups.push_back(gr);
            // name and version
            _basecall_group_descriptions[gr] = detect_basecall_group_id(gr);
            auto & bc_desc = _basecall_group_descriptions.at(gr);
            // subgroups
            for (unsigned st = 0; st < 3; ++st)
            {
                bc_desc.have_subgroup[st] =
                    Base::group_exists(basecall_strand_group_path(gr, st));
                if (bc_desc.have_subgroup[st])
                {
                    _basecall_strand_groups[st].push_back(gr);
                    // fastq
                    bc_desc.have_fastq[st] =
                        have_basecall_fastq_unpack(st, gr) or
                        have_basecall_fastq_pack(st, gr);
                    // events
                    bc_desc.have_events[st] =
                        have_basecall_events_unpack(st, gr) or
                        have_basecall_events_pack(st, gr);
                    if (st == 0)
                    {
                        // ed_gr
                        bc_desc.ed_gr = detect_basecall_eventdetection_group(gr);
                    }
                    if (st == 2)
                    {
                        // alignment
                        bc_desc.have_alignment =
                            have_basecall_alignment_unpack(gr)
                            or have_basecall_alignment_pack(gr);
                    }
                }
            }
            // bc_1d_gr
            if (bc_desc.have_subgroup[0] or bc_desc.have_subgroup[1])
            {
                bc_desc.bc_1d_gr = gr;
            }
            else if (bc_desc.have_subgroup[2])
            {
                bc_desc.bc_1d_gr = detect_basecall_1d_group(gr);
            }
            // model
            for (unsigned st = 0; st < 2; ++st)
            {
                bc_desc.have_model[st] =
                    not bc_desc.bc_1d_gr.empty()
                    and Base::dataset_exists(basecall_model_path(bc_desc.bc_1d_gr, st));
            }
        }
    }
    Basecall_Group_Description
    detect_basecall_group_id(std::string const & gr) const
    {
        Basecall_Group_Description res;
        res.name = "?";
        res.version = "?";
        auto am = get_basecall_params(gr);
        if (am.count("name"))
        {
            if (am.at("name") == "ONT Sequencing Workflow")
            {
                res.name = "metrichor";
                res.version = (am.count("chimaera version")? am.at("chimaera version") : "?") + "+" +
                    (am.count("dragonet version")? am.at("dragonet version") : "?");
            }
            else if (am.at("name") == "MinKNOW-Live-Basecalling")
            {
                res.name = "minknow";
                res.version = (am.count("version")? am.at("version") : "?");
            }
            else if (am.at("name") == "ONT Albacore Sequencing Software")
            {
                res.name = "albacore";
                res.version = (am.count("version")? am.at("version") : "?");
            }
        }
        return res;
    }
    std::string
    detect_basecall_1d_group(std::string const & gr) const
    {
        std::string path = basecall_group_path(gr) + "/basecall_1d";
        if (Base::attribute_exists(path))
        {
            std::string tmp;
            Base::read(path, tmp);

            // Metrichor writes "Analyses" before the group, Albacore does not
            std::string possible_prefix = "Analyses";
            std::string pref = "";

            if(tmp.substr(0, possible_prefix.length()) == possible_prefix) {
                // metrichor version
                pref = possible_prefix + "/" + basecall_group_prefix();
            } else {
                // albacore version
                pref = basecall_group_prefix();
            }

            if (tmp.size() >= pref.size()
                and tmp.substr(0, pref.size()) == pref)
            {
                auto gr_1d = tmp.substr(pref.size());
                if (have_basecall_group(gr_1d))
                {
                    return gr_1d;
                }
            }
        }
        return gr;
    }
    std::string
    detect_basecall_eventdetection_group(std::string const & gr) const
    {
        auto bc_params = get_basecall_params(gr);
        if (bc_params.count("event_detection"))
        {
            auto && tmp = bc_params.at("event_detection");
            auto pref = eventdetection_root_path().substr(1) + "/" + eventdetection_group_prefix();
            if (tmp.substr(0, pref.size()) == pref)
            {
                auto ed_gr = tmp.substr(pref.size());
                if (have_eventdetection_group(ed_gr))
                {
                    return ed_gr;
                }
            }
        }
        if (have_basecall_events_pack(0, gr))
        {
            auto ev_pack = get_basecall_events_pack(0, gr);
            auto ed_gr = ev_pack.ed_gr;
            if (have_eventdetection_group(ed_gr))
            {
                return ed_gr;
            }
        }
        return "";
    }
    double
    get_basecall_median_sd_temp(std::string const & gr) const
    {
        std::string segmentation_link_path = basecall_group_path(gr) + "/segmentation";
        if (not Base::attribute_exists(segmentation_link_path)) return 0.0;
        std::string segmentation_path;
        Base::read(segmentation_link_path, segmentation_path);
        std::string median_sd_temp_path = "/" + segmentation_path + "/Summary/split_hairpin/median_sd_temp";
        if (not Base::attribute_exists(median_sd_temp_path)) return 0.0;
        double res;
        Base::read(median_sd_temp_path, res);
        return res;
    }

    //
    // Functions that fill in empty arguments with default values
    //
    std::string const &
    fill_raw_samples_read_name(std::string const & rn) const
    {
        return (not rn.empty() or _raw_samples_read_names.empty()
                ? rn
                : _raw_samples_read_names.front());
    }
    std::string const &
    fill_eventdetection_group(std::string const & gr) const
    {
        return (not gr.empty() or _eventdetection_groups.empty()
                ? gr
                : _eventdetection_groups.front());
    }
    std::string const &
    fill_eventdetection_read_name(std::string const & gr, std::string const & rn) const
    {
        return (not rn.empty()
                or _eventdetection_read_names.count(gr) == 0
                or _eventdetection_read_names.at(gr).empty()
                ? rn
                : _eventdetection_read_names.at(gr).front());
    }
    std::string const &
    fill_basecall_group(unsigned st, std::string const & gr) const
    {
        return (not gr.empty()
                or _basecall_strand_groups.at(st).empty()
                ? gr
                : _basecall_strand_groups.at(st).front());
    }
    std::string const &
    fill_basecall_1d_group(unsigned st, std::string const & gr) const
    {
        auto && _gr = fill_basecall_group(st, gr);
        return get_basecall_1d_group(_gr);
    }

    //
    // Packing interface
    //
    Raw_Samples_Pack
    get_raw_samples_pack(std::string const & rn) const
    {
        Raw_Samples_Pack rs_pack;
        auto path = raw_samples_pack_path(rn);
        rs_pack.read(*this, path);
        return rs_pack;
    }
    void
    add_raw_samples(std::string const & rn, Raw_Samples_Pack const & rs_pack)
    {
        auto path = raw_samples_pack_path(rn);
        rs_pack.write(*this, path);
        reload();
    }
    Raw_Int_Samples_Dataset
    get_raw_int_samples_dataset(std::string const & rn = std::string()) const
    {
        Raw_Int_Samples_Dataset res;
        auto && _rn = fill_raw_samples_read_name(rn);
        res.first = get_raw_int_samples(_rn);
        res.second = get_raw_samples_params(_rn);
        return res;
    }
    Raw_Samples_Dataset
    get_raw_samples_dataset(std::string const & rn = std::string()) const
    {
        Raw_Samples_Dataset res;
        auto && _rn = fill_raw_samples_read_name(rn);
        res.first = get_raw_samples(_rn);
        res.second = get_raw_samples_params(_rn);
        return res;
    }
    void
    add_raw_samples_dataset(std::string const & rn, Raw_Int_Samples_Dataset const & rsi_ds)
    {
        add_raw_samples(rn, rsi_ds.first);
        add_raw_samples_params(rn, rsi_ds.second);
    }
    EventDetection_Events_Pack
    get_eventdetection_events_pack(
        std::string const & gr, std::string const & rn) const
    {
        EventDetection_Events_Pack ede_pack;
        ede_pack.read(*this, eventdetection_events_pack_path(gr, rn));
        return ede_pack;
    }
    void
    add_eventdetection_events(
        std::string const & gr, std::string const & rn,
        EventDetection_Events_Pack const & ede_pack)
    {
        ede_pack.write(*this, eventdetection_events_pack_path(gr, rn));
        reload();
    }
    EventDetection_Events_Dataset
    get_eventdetection_events_dataset(
        std::string const & gr, std::string const & rn) const
    {
        EventDetection_Events_Dataset ede_ds;
        ede_ds.first = get_eventdetection_events(gr, rn);
        ede_ds.second = get_eventdetection_events_params(gr, rn);
        return ede_ds;
    }
    void
    add_eventdetection_events_dataset(
        std::string const & gr, std::string const & rn,
        EventDetection_Events_Dataset const & ede_ds)
    {
        add_eventdetection_events(gr, rn, ede_ds.first);
        add_eventdetection_events_params(gr, rn, ede_ds.second);
    }
    //
    Basecall_Fastq_Pack
    get_basecall_fastq_pack(unsigned st, std::string const & gr) const
    {
        Basecall_Fastq_Pack fq_pack;
        auto p = basecall_fastq_pack_path(gr, st);
        fq_pack.read(*this, p);
        return fq_pack;        
    }
    void
    add_basecall_fastq(unsigned st, std::string const & gr, Basecall_Fastq_Pack const & fq_pack)
    {
        auto p = basecall_fastq_pack_path(gr, st);
        fq_pack.write(*this, p);
        reload();
    }
    //
    Basecall_Events_Pack
    get_basecall_events_pack(unsigned st, std::string const & gr) const
    {
        auto p = basecall_events_pack_path(gr, st);
        Basecall_Events_Pack ev_pack;
        ev_pack.read(*this, p);
        return ev_pack;
    }
    void
    add_basecall_events(unsigned st, std::string const & gr, Basecall_Events_Pack const & ev_pack)
    {
        auto p = basecall_events_pack_path(gr, st);
        ev_pack.write(*this, p);
        reload();
    }
    Basecall_Events_Dataset
    get_basecall_events_dataset(unsigned st, std::string const & gr) const
    {
        Basecall_Events_Dataset bce_ds;
        bce_ds.first = get_basecall_events(st, gr);
        bce_ds.second = get_basecall_events_params(st, gr);
        return bce_ds;
    }
    void
    add_basecall_events_dataset(unsigned st, std::string const & gr, Basecall_Events_Dataset const & bce_ds)
    {
        add_basecall_events(st, gr, bce_ds.first);
        add_basecall_events_params(st, gr, bce_ds.second);
    }
    //
    Basecall_Alignment_Pack
    get_basecall_alignment_pack(std::string const & gr) const
    {
        Basecall_Alignment_Pack al_pack;
        auto p = basecall_alignment_pack_path(gr);
        al_pack.read(*this, p);
        return al_pack;
    }
    void
    add_basecall_alignment(std::string const & gr, Basecall_Alignment_Pack const & al_pack)
    {
        auto p = basecall_alignment_pack_path(gr);
        al_pack.write(*this, p);
        reload();
    }

    //
    // Packers & Unpackers
    //
    static Raw_Samples_Pack
    pack_rw(Raw_Int_Samples_Dataset const & rsi_ds)
    {
        Raw_Samples_Pack rsp;
        rsp.params = rsi_ds.second;
        std::tie(rsp.signal, rsp.signal_params) = rw_coder().encode(rsi_ds.first, true);
        return rsp;
    }
    static Raw_Int_Samples_Dataset
    unpack_rw(Raw_Samples_Pack const & rs_pack)
    {
        Raw_Int_Samples_Dataset rsi_ds;
        rsi_ds.second = rs_pack.params;
        rsi_ds.first = rw_coder().decode< Raw_Int_Sample >(rs_pack.signal, rs_pack.signal_params);
        return rsi_ds;
    }
    static std::pair< std::vector< long long >, std::vector< long long > >
    pack_event_start_length(
        unsigned num_events,
        std::function< long long(unsigned) > get_start,
        std::function< long long(unsigned) > get_length,
        long long start_time)
    {
        std::pair< std::vector< long long >, std::vector< long long > > res;
        auto & skip = res.first;
        auto & len = res.second;
        for (unsigned i = 0; i < num_events; ++i)
        {
            auto si = get_start(i);
            auto li = get_length(i);
            skip.push_back(si - start_time);
            len.push_back(li);
            start_time = si + li;
        }
        return res;
    }
    static void
    unpack_event_start_length(
        std::vector< long long > const & skip,
        std::vector< long long > const & len,
        std::function< void(unsigned, long long) > set_start,
        std::function< void(unsigned, long long) > set_length,
        long long start_time)
    {
        for (unsigned i = 0; i < skip.size(); ++i)
        {
            auto si = start_time + skip[i];
            auto li = len[i];
            set_start(i, si);
            set_length(i, li);
            start_time = si + li;
        }
    }
    static void
    unpack_event_mean_stdv(
        unsigned num_events,
        std::function< long long(unsigned) > get_start,
        std::function< long long(unsigned) > get_length,
        std::function< void(unsigned, double) > set_mean,
        std::function< void(unsigned, double) > set_stdv,
        std::vector< Raw_Sample > const & rs,
        long long rs_start_time,
        int offset)
    {
        for (unsigned i = 0; i < num_events; ++i)
        {
            long long rs_start_idx = get_start(i) - rs_start_time + offset;
            long long rs_end_idx = rs_start_idx + get_length(i);
            if (i == 0 and rs_start_idx < 0) rs_start_idx = 0;
            if (i == num_events - 1 and rs_end_idx > (long long)rs.size()) rs_end_idx = rs.size();
            if (rs_start_idx < 0
                or rs_end_idx <= rs_start_idx
                or rs_end_idx > (long long)rs.size())
            {
                LOG_THROW
                    << "bad index: rs_start_idx=" << rs_start_idx
                    << " rs_end_idx=" << rs_end_idx
                    << " i=" << i
                    << " length(i)=" << get_length(i)
                    << " rs_size=" << rs.size()
                    << " offset=" << offset;
            }
            bool all_equal = true;
            double s = 0.0;
            double s2 = 0.0;
            unsigned n = rs_end_idx - rs_start_idx;
            for (unsigned j = 0; j < n; ++j)
            {
                double x = rs[rs_start_idx + j];
                if (j > 0 and all_equal)
                {
                    all_equal = rs[rs_start_idx + j] == rs[rs_start_idx];
                }
                s += x;
                s2 += x * x;
            }
            set_mean(i, s / n);
            if (n > 1 and not all_equal)
            {
                double x = (s2 - s*s/n)/n;
                set_stdv(i, x > 1e-3? std::sqrt(x) : 0);
            }
            else
            {
                set_stdv(i, 0);
            }
        }
    }
    static EventDetection_Events_Pack
    pack_ed(EventDetection_Events_Dataset const & ede_ds)
    {
        EventDetection_Events_Pack ede_pack;
        auto & ede = ede_ds.first;
        auto & ede_params = ede_ds.second;
        ede_pack.params = ede_params;
        std::vector< long long > skip;
        std::vector< long long > len;
        std::tie(skip, len) = pack_event_start_length(
            ede.size(),
            [&] (unsigned i) { return ede.at(i).start; },
            [&] (unsigned i) { return ede.at(i).length; },
            ede_params.start_time);
        std::tie(ede_pack.skip, ede_pack.skip_params) = ed_skip_coder().encode(skip, false);
        std::tie(ede_pack.len, ede_pack.len_params) = ed_len_coder().encode(len, false);
        return ede_pack;
    }
    static EventDetection_Events_Dataset
    unpack_ed(EventDetection_Events_Pack const & ede_pack,
              Raw_Samples_Dataset const & rs_ds)
    {
        EventDetection_Events_Dataset res;
        auto & ede_params = ede_pack.params;
        auto & rs = rs_ds.first;
        auto & rs_params = rs_ds.second;
        res.second = ede_params;
        auto skip = ed_skip_coder().decode< long long >(ede_pack.skip, ede_pack.skip_params);
        auto len = ed_len_coder().decode< long long >(ede_pack.len, ede_pack.len_params);
        if (skip.size() != len.size())
        {
            LOG_THROW
                << "wrong dataset size: skip_size=" << skip.size()
                << " len_size=" << len.size();
        }
        auto & ede = res.first;
        ede.resize(skip.size());
        unpack_event_start_length(
            skip,
            len,
            [&] (unsigned i, long long x) { return ede.at(i).start = x; },
            [&] (unsigned i, long long x) { return ede.at(i).length = x; },
            ede_params.start_time);
        int offset = 0;
        static bool warned = false;
        if (offset != 0 and not warned)
        {
            LOG(warning) << "using workaround for old off-by-one ed events bug\n";
            warned = true;
        }
        unpack_event_mean_stdv(
            ede.size(),
            [&] (unsigned i) { return ede.at(i).start; },
            [&] (unsigned i) { return ede.at(i).length; },
            [&] (unsigned i, double x) { return ede.at(i).mean = x; },
            [&] (unsigned i, double x) { return ede.at(i).stdv = x; },
            rs,
            rs_params.start_time,
            offset);
        return res;
    }
    static Basecall_Fastq_Pack
    pack_fq(std::string const & fq, unsigned qv_bits = 5)
    {
        static unsigned const max_qv_bits = 5;
        static std::uint8_t const max_qv = ((std::uint8_t)1 << max_qv_bits) - 1;
        Basecall_Fastq_Pack fq_pack;
        auto fqa = split_fq(fq);
        fq_pack.read_name = fqa[0];
        std::vector< std::int8_t > bp(fqa[1].begin(), fqa[1].end());
        qv_bits = std::min(qv_bits, max_qv_bits);
        auto qv_mask = max_qv & (max_qv << (max_qv_bits - qv_bits));
        fq_pack.qv_bits = qv_bits;
        std::vector< std::uint8_t > qv;
        for (auto c : fqa[3])
        {
            std::uint8_t val = (std::uint8_t)(c - 33);
            val = std::min(val, max_qv);
            val &= qv_mask;
            qv.push_back(val);
        }
        std::tie(fq_pack.bp, fq_pack.bp_params) = fq_bp_coder().encode(bp, false);
        std::tie(fq_pack.qv, fq_pack.qv_params) = fq_qv_coder().encode(qv, false);
        return fq_pack;
    }
    static std::string
    unpack_fq(Basecall_Fastq_Pack const & fq_pack)
    {
        std::string res;
        res += "@";
        res += fq_pack.read_name;
        res += "\n";
        auto bp = fq_bp_coder().decode< std::int8_t >(fq_pack.bp, fq_pack.bp_params);
        for (auto c : bp) res += c;
        res += "\n+\n";
        auto qv = fq_qv_coder().decode< std::uint8_t >(fq_pack.qv, fq_pack.qv_params);
        for (auto c : qv) res += (char)33 + c;
        res += "\n";
        return res;
    }
    static Basecall_Events_Pack
    pack_ev(Basecall_Events_Dataset const & ev_ds,
            Basecall_Group_Description const & bc_desc,
            std::string const & sq,
            std::vector< EventDetection_Event > const & ed,
            std::string const & ed_gr,
            Channel_Id_Params const & cid_params,
            double median_sd_temp,
            unsigned p_model_state_bits)
    {
        Basecall_Events_Pack ev_pack;
        ev_pack.params = ev_ds.second;
        auto & ev = ev_ds.first;
        ev_pack.name = bc_desc.name;
        ev_pack.version = bc_desc.version;
        ev_pack.ed_gr = ed_gr;
        ev_pack.start_time = time_to_int(ev[0].start, cid_params);
        ev_pack.state_size = ev[0].get_model_state().size();
        ev_pack.median_sd_temp = median_sd_temp;
        ev_pack.p_model_state_bits = p_model_state_bits;
        std::vector< long long > rel_skip;
        std::vector< long long > skip;
        std::vector< long long > len;
        std::vector< std::uint8_t > mv;
        std::vector< std::uint16_t > p_model_state;
        // first pack start/duration
        if (not ed_gr.empty())
        {
            // pack relative to ed events
            long long j = -1;
            for (unsigned i = 0; i < ev.size(); ++i)
            {
                auto ti = time_to_int(ev[i].start, cid_params);
                auto last_j = j++;
                while (j < (long long)ed.size() and ed[j].start < ti) ++j;
                if (j == (long long)ed.size())
                {
                    LOG_THROW
                        << "no matching ed event: i=" << i
                        << " ev[i]=(" << ti
                        << "," << time_to_int(ev[i].length, cid_params)
                        << "," << ev[i].mean
                        << "," << ev[i].stdv
                        << ")";
                }
                rel_skip.push_back(j - last_j - 1);
            }
            std::tie(ev_pack.rel_skip, ev_pack.rel_skip_params) = ev_rel_skip_coder().encode(rel_skip, false);
        }
        else
        {
            // pack start&length as for ed events
            std::tie(skip, len) = pack_event_start_length(
                ev.size(),
                [&] (unsigned i) { return time_to_int(ev.at(i).start, cid_params); },
                [&] (unsigned i) { return time_to_int(ev.at(i).length, cid_params); },
                ev_pack.start_time);
            std::tie(ev_pack.skip, ev_pack.skip_params) = ed_skip_coder().encode(skip, false);
            std::tie(ev_pack.len, ev_pack.len_params) = ed_len_coder().encode(len, false);
        }
        unsigned sq_pos = 0;
        for (unsigned i = 0; i < ev.size(); ++i)
        {
            auto s = ev[i].get_model_state();
            if (s.size() != ev_pack.state_size)
            {
                LOG_THROW
                    << "unexpected state size: i=" << i
                    << " s=" << s
                    << " expected_size=" << ev_pack.state_size;
            }
            // check if move is valid
            if (ev[i].move < 0 or ev[i].move > std::numeric_limits< uint8_t >::max())
            {
                LOG_THROW
                    << "invalid move: i=" << i
                    << "ev[i].move=" << ev[i].move;
            }
            int real_move = ev[i].move;
            if (sq.substr(sq_pos + real_move, ev_pack.state_size) != s)
            {
                // move is not valid, compute alternative:
                // allow move > state_size only if previous state is homopolymer
                auto next_sq_pos = sq.find(s, sq_pos);
                if (next_sq_pos != std::string::npos
                    and (next_sq_pos <= sq_pos + ev_pack.state_size
                         or sq.substr(sq_pos, ev_pack.state_size) == std::string(ev_pack.state_size, sq[sq_pos])))
                {
                    real_move = next_sq_pos - sq_pos;
                }
                else
                {
                    real_move = -1;
                }
                if (real_move >= 0)
                {
                    LOG(warning)
                        << "using workaround for invalid move: i=" << i
                        << " sq=" << sq.substr(sq_pos, 2 * ev_pack.state_size)
                        << " move[i]=" << ev[i].move
                        << " state[i]=" << s
                        << " real_move=" << real_move << std::endl;
                }
                else
                {
                    LOG_THROW
                        << "invalid move: i=" << i
                        << " sq=" << sq.substr(sq_pos, 2 * ev_pack.state_size)
                        << " move[i]=" << ev[i].move
                        << " state[i]=" << s;
                }
            }
            mv.push_back(real_move);
            sq_pos += real_move;
            // p_model_state
            std::uint16_t p_model_state_val = ev[i].p_model_state * (1u << p_model_state_bits);
            if (p_model_state_val >= (1u << p_model_state_bits)) p_model_state_val = (1u << p_model_state_bits) - 1;
            p_model_state.push_back(p_model_state_val);
        }
        if (sq_pos + ev_pack.state_size != sq.size())
        {
            LOG_THROW
                << "leftover base sequence: sq_size=" << sq.size()
                << " sq_end_pos=" << sq_pos + ev_pack.state_size;
        }
        std::tie(ev_pack.move, ev_pack.move_params) = ev_move_coder().encode(mv, false);
        std::tie(ev_pack.p_model_state, ev_pack.p_model_state_params) = bit_packer().encode(p_model_state, p_model_state_bits);
        return ev_pack;
    } // pack_ev()
    static std::vector< EventDetection_Event >
    unpack_implicit_ed(Basecall_Events_Pack const & ev_pack,
                       Raw_Samples_Dataset const & rs_ds)
    {
        std::vector< EventDetection_Event > ede;
        auto & rs = rs_ds.first;
        auto & rs_params = rs_ds.second;
        auto skip = ed_skip_coder().decode< long long >(ev_pack.skip, ev_pack.skip_params);
        auto len = ed_len_coder().decode< long long >(ev_pack.len, ev_pack.len_params);
        if (skip.empty() or skip.size() != len.size())
        {
            LOG_THROW
                << "wrong dataset size: skip_size=" << skip.size()
                << " len_size=" << len.size();
        }
        ede.resize(skip.size());
        unpack_event_start_length(
            skip,
            len,
            [&] (unsigned i, long long x) { return ede.at(i).start = x; },
            [&] (unsigned i, long long x) { return ede.at(i).length = x; },
            ev_pack.start_time);
        int offset = 0;
        static bool warned = false;
        if (offset != 0 and not warned)
        {
            LOG(warning) << "using workaround for bug in "
                         << ev_pack.name << ":" << ev_pack.version << "\n";
            warned = true;
        }
        unpack_event_mean_stdv(
            ede.size(),
            [&] (unsigned i) { return ede.at(i).start; },
            [&] (unsigned i) { return ede.at(i).length; },
            [&] (unsigned i, double x) { return ede.at(i).mean = x; },
            [&] (unsigned i, double x) { return ede.at(i).stdv = x; },
            rs,
            rs_params.start_time,
            offset);
        return ede;
    }
    static Basecall_Events_Dataset
    unpack_ev(Basecall_Events_Pack const & ev_pack,
              std::string const & sq,
              std::vector< EventDetection_Event > const & ed,
              Channel_Id_Params const & cid_params)
    {
        Basecall_Events_Dataset ev_ds;
        ev_ds.second = ev_pack.params;
        auto & ev = ev_ds.first;
        std::vector< long long > rel_skip;
        if (not ev_pack.rel_skip.empty())
        {
            rel_skip = ev_rel_skip_coder().decode< long long >(ev_pack.rel_skip, ev_pack.rel_skip_params);
        }
        auto mv = ev_move_coder().decode< std::uint8_t >(ev_pack.move, ev_pack.move_params);
        auto p_model_state = bit_packer().decode< std::uint16_t >(ev_pack.p_model_state, ev_pack.p_model_state_params);
        if ((not rel_skip.empty() and rel_skip.size() != mv.size()) or p_model_state.size() != mv.size())
        {
            LOG_THROW
                << "wrong dataset size: rel_skip_size=" << rel_skip.size()
                << " mv_size=" << mv.size()
                << " p_model_state_size=" << p_model_state.size();
        }
        ev.resize(mv.size());
        long long j = -1;
        std::string s;
        unsigned sq_pos = 0;
        unsigned p_model_state_bits;
        std::istringstream(ev_pack.p_model_state_params.at("num_bits")) >> p_model_state_bits;
        long long unsigned max_p_model_state_int = 1llu << p_model_state_bits;
        for (unsigned i = 0; i < ev.size(); ++i)
        {
            j += (not rel_skip.empty()? rel_skip[i] : 0) + 1;
            ev[i].start = time_to_float(ed[j].start, cid_params);
            ev[i].length = time_to_float(ed[j].length, cid_params);
            ev[i].mean = ed[j].mean;
            ev[i].stdv = ed[j].stdv;
            if (ev[i].stdv == 0.0) ev[i].stdv = ev_pack.median_sd_temp;
            ev[i].move = mv[i];
            if (i > 0) s = s.substr(mv[i]); // apply move
            while (s.size() < ev_pack.state_size) s += sq[sq_pos++];
            std::copy(s.begin(), s.end(), ev[i].model_state.begin());
            if (ev_pack.state_size < MAX_K_LEN) ev[i].model_state[ev_pack.state_size] = 0;
            ev[i].p_model_state = (double)p_model_state[i] / max_p_model_state_int;
        }
        return ev_ds;
    } // unpack_ev()
    static Basecall_Alignment_Pack
    pack_al(std::vector< Basecall_Alignment_Entry > const & al,
            std::string const & sq)
    {
        Basecall_Alignment_Pack al_pack;
        std::array< std::vector< uint8_t > , 2 > step_v;
        std::vector< int8_t > mv;
        step_v[0].reserve(al.size());
        step_v[1].reserve(al.size());
        mv.reserve(al.size());
        std::array< int, 2 > start_index = {{ -1, -1 }};
        std::array< int, 2 > next_index = {{ -1, -1 }};
        std::array< int, 2 > delta = {{ 1, -1 }};
        auto get_idx = [&] (unsigned i, unsigned k) {
            return k == 0? al[i].template_index : al[i].complement_index;
        };
        unsigned pos = 0;
        for (unsigned i = 0; i < al.size(); ++i)
        {
            for (unsigned k = 0; k < 2; ++k)
            {
                auto idx = get_idx(i, k);
                if (idx >= 0)
                {
                    if (start_index[k] < 0)
                    {
                        start_index[k] = idx;
                        next_index[k] = idx;
                    }
                    if (idx != next_index[k])
                    {
                        LOG_THROW
                            << "bad index: idx=" << idx
                            << " next_index=" << next_index[k];
                    }
                    step_v[k].push_back(1);
                    next_index[k] += delta[k];
                }
                else // idx < 0
                {
                    step_v[k].push_back(0);
                }
            }
            // compute move
            auto kmer = al[i].get_kmer();
            size_t next_pos = sq.find(kmer, pos);
            if (next_pos == std::string::npos)
            {
                LOG_THROW
                    << "missing kmer in 2d seq";
            }
            if (next_pos - pos > std::numeric_limits< int8_t >::max())
            {
                LOG_THROW
                    << "bad move: next_pos=" << next_pos
                    << " pos=" << pos;
            }
            mv.push_back(next_pos - pos);
            pos = next_pos;
        }
        if (start_index[0] < 0)
        {
            LOG_THROW
                << "no template events";
        }
        if (start_index[1] < 0)
        {
            LOG_THROW
                << "no complement events";
        }
        al_pack.template_index_start = start_index[0];
        al_pack.complement_index_start = start_index[1];
        al_pack.kmer_size = al[0].get_kmer().size();
        std::tie(al_pack.template_step, al_pack.template_step_params) = bit_packer().encode(step_v[0], 1);
        std::tie(al_pack.complement_step, al_pack.complement_step_params) = bit_packer().encode(step_v[1], 1);
        std::tie(al_pack.move, al_pack.move_params) = ev_move_coder().encode(mv, false);
        return al_pack;
    } // pack_al()
    static std::vector< Basecall_Alignment_Entry >
    unpack_al(Basecall_Alignment_Pack const & al_pack,
              std::string const & sq)
    {
        std::vector< Basecall_Alignment_Entry > al;
        std::array< std::vector< uint8_t >, 2 > step_v =
            {{ bit_packer().decode< uint8_t >(al_pack.template_step, al_pack.template_step_params),
               bit_packer().decode< uint8_t >(al_pack.complement_step, al_pack.complement_step_params) }};
        auto mv = ev_move_coder().decode< int8_t >(al_pack.move, al_pack.move_params);
        if (step_v[1].size() != step_v[0].size()
            or mv.size() != step_v[0].size())
        {
            LOG_THROW
                << "wrong dataset size: step_v[0]_size=" << step_v[0].size()
                << " step_v[1]_size=" << step_v[1].size()
                << " mv_size=" << mv.size();
        }
        al.resize(step_v[0].size());
        std::array< unsigned, 2 > crt_index = {{ al_pack.template_index_start, al_pack.complement_index_start }};
        std::array< int, 2 > delta = {{ 1, -1 }};
        auto pos = 0;
        auto set_idx = [&] (unsigned i, unsigned k, int val) {
            if (k == 0)
            {
                al[i].template_index = val;
            }
            else
            {
                al[i].complement_index = val;
            }
        };
        for (unsigned i = 0; i < step_v[0].size(); ++i)
        {
            for (unsigned k = 0; k < 2; ++k)
            {
                if (step_v[k][i] > 0)
                {
                    set_idx(i, k, crt_index[k]);
                    crt_index[k] += delta[k];
                }
                else
                {
                    set_idx(i, k, -1);
                }
            }
            // set kmer
            pos += mv[i];
            std::copy(sq.begin() + pos, sq.begin() + pos + al_pack.kmer_size, al[i].kmer.begin());
            if (al_pack.kmer_size < MAX_K_LEN) al[i].kmer[al_pack.kmer_size] = 0;
        }
        return al;
    } // unpack_al()

    //
    // Fast5 internal paths
    //
    static std::string file_version_path() { return "/file_version"; }
    static std::string channel_id_path()   { return "/UniqueGlobalKey/channel_id"; }
    static std::string tracking_id_path()  { return "/UniqueGlobalKey/tracking_id"; }
    static std::string context_tags_path()  { return "/UniqueGlobalKey/context_tags"; }
    static std::string sequences_path()    { return "/Sequences/Meta"; }
    static std::string raw_samples_root_path() { return "/Raw/Reads"; }
    static std::string raw_samples_params_path(std::string const & rn)
    {
        return raw_samples_root_path() + "/" + rn;
    }
    static std::string raw_samples_path(std::string const & rn)
    {
        return raw_samples_root_path() + "/" + rn + "/Signal";
    }
    static std::string raw_samples_pack_path(std::string const & rn)
    {
        return raw_samples_path(rn) + "_Pack";
    }
    static std::string raw_samples_params_pack_path(std::string const & rn)
    {
        return raw_samples_pack_path(rn) + "/params";
    }
    static std::string eventdetection_root_path() { return "/Analyses"; }
    static std::string eventdetection_group_prefix() { return "EventDetection_"; }
    static std::string eventdetection_group_path(std::string const & gr)
    {
        return eventdetection_root_path() + "/" + eventdetection_group_prefix() + gr;
    }
    static std::string eventdetection_events_params_path(std::string const & gr, std::string const & rn)
    {
        return eventdetection_group_path(gr) + "/Reads/" + rn;
    }
    static std::string eventdetection_events_path(std::string const & gr, std::string const & rn)
    {
        return eventdetection_group_path(gr) + "/Reads/" + rn + "/Events";
    }
    static std::string eventdetection_events_pack_path(std::string const & gr, std::string const & rn)
    {
        return eventdetection_events_path(gr, rn) + "_Pack";
    }
    static std::string eventdetection_events_params_pack_path(std::string const & gr, std::string const & rn)
    {
        return eventdetection_events_pack_path(gr, rn) + "/params";
    }
    static std::string basecall_root_path() { return "/Analyses"; }
    static std::string basecall_group_prefix() { return "Basecall_"; }
    static std::string strand_name(unsigned st)
    {
        static const std::array< std::string, 3 > _strand_name =
            {{ "template", "complement", "2D" }};
        return _strand_name.at(st);
    }
    static std::string basecall_strand_subgroup(unsigned st)
    {
        return std::string("BaseCalled_") + strand_name(st);
    }
    static std::string basecall_group_path(std::string const & gr)
    {
        return basecall_root_path() + "/" + basecall_group_prefix() + gr;
    }
    static std::string basecall_strand_group_path(std::string const & gr, unsigned st)
    {
        return basecall_group_path(gr) + "/" + basecall_strand_subgroup(st);
    }
    static std::string basecall_log_path(std::string const & gr)
    {
        return basecall_group_path(gr) + "/Log";
    }
    static std::string basecall_fastq_path(std::string const & gr, unsigned st)
    {
        return basecall_strand_group_path(gr, st) + "/Fastq";
    }
    static std::string basecall_fastq_pack_path(std::string const & gr, unsigned st)
    {
        return basecall_fastq_path(gr, st) + "_Pack";
    }
    static std::string basecall_model_path(std::string const & gr, unsigned st)
    {
        return basecall_strand_group_path(gr, st) + "/Model";
    }
    static std::string basecall_model_file_path(std::string const & gr, unsigned st)
    {
        return basecall_group_path(gr) + "/Summary/basecall_1d_" + strand_name(st) + "/model_file";
    }
    static std::string basecall_events_path(std::string const & gr, unsigned st)
    {
        return basecall_strand_group_path(gr, st) + "/Events";
    }
    static std::string basecall_events_pack_path(std::string const & gr, unsigned st)
    {
        return basecall_events_path(gr, st) + "_Pack";
    }
    static std::string basecall_events_params_pack_path(std::string const & gr, unsigned st)
    {
        return basecall_events_pack_path(gr, st) + "/params";
    }
    static std::string basecall_alignment_path(std::string const & gr)
    {
        return basecall_strand_group_path(gr, 2) + "/Alignment";
    }
    static std::string basecall_alignment_pack_path(std::string const & gr)
    {
        return basecall_alignment_path(gr) + "_Pack";
    }
    static std::string basecall_config_path(std::string const & gr)
    {
        return basecall_group_path(gr) + "/Configuration";
    }
    static std::string basecall_summary_path(std::string const & gr)
    {
        return basecall_group_path(gr) + "/Summary";
    }
    //
    // Packers
    //
    static Huffman_Packer const & rw_coder()          { return Huffman_Packer::get_coder("fast5_rw_1"); }
    static Huffman_Packer const & ed_skip_coder()     { return Huffman_Packer::get_coder("fast5_ed_skip_1"); }
    static Huffman_Packer const & ed_len_coder()      { return Huffman_Packer::get_coder("fast5_ed_len_1"); }
    static Huffman_Packer const & fq_bp_coder()       { return Huffman_Packer::get_coder("fast5_fq_bp_1"); }
    static Huffman_Packer const & fq_qv_coder()       { return Huffman_Packer::get_coder("fast5_fq_qv_1"); }
    static Huffman_Packer const & ev_rel_skip_coder() { return Huffman_Packer::get_coder("fast5_ev_rel_skip_1"); }
    static Huffman_Packer const & ev_move_coder()     { return Huffman_Packer::get_coder("fast5_ev_move_1"); }
    static Bit_Packer     const & bit_packer()        { return Bit_Packer::get_packer(); }
}; // class File

} // namespace fast5

#endif
