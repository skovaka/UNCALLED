#
# Part of: https://github.com/mateidavid/fast5
#
# (c) 2017: Matei David, Ontario Institute for Cancer Research
# MIT License
#

from cython.operator cimport dereference as deref

from libc.stdint cimport int16_t
from libcpp cimport bool
from libcpp.map cimport map as cmap
from libcpp.memory cimport unique_ptr
from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "fast5.hpp" namespace "fast5":

    cdef string cpp_version "fast5::version"

    cppclass Cpp_Logger "logger::Logger":
        @staticmethod
        void set_levels_from_options(vector[string]) except +

    ctypedef cmap[string, string] Attr_Map

    struct Channel_Id_Params:
        string channel_number
        double digitisation
        double offset
        double range
        double sampling_rate

    ctypedef Attr_Map Tracking_Id_Params

    ctypedef Attr_Map Sequences_Params

    struct Raw_Samples_Params:
        string read_id
        long long read_number
        long long start_mux
        long long start_time
        long long duration

    ctypedef float Raw_Sample

    ctypedef int16_t Raw_Int_Sample

    struct EventDetection_Events_Params:
        string read_id
        long long read_number
        long long scaling_used
        long long start_mux
        long long start_time
        long long duration
        double median_before
        unsigned abasic_found

    struct EventDetection_Event:
        double mean
        double stdv
        long long start
        long long length

    struct Basecall_Model_Params:
        double scale
        double shift
        double drift
        double var
        double scale_sd
        double var_sd

    struct Basecall_Model_State:
        double level_mean
        double level_stdv
        double sd_mean
        double sd_stdv
        #char kmer[8]

    struct Basecall_Events_Params:
        double start_time
        double duration

    struct Basecall_Event:
        double mean
        double stdv
        double start
        double length
        double p_model_state
        long long move
        #char model_state[8]

    struct Basecall_Alignment_Entry:
        long long template_index
        long long complement_index
        #char kmer[8]

    struct Basecall_Group_Description:
        string name
        string version
        string ed_gr
        string bc_1d_gr
        bool have_subgroup[3]
        bool have_fastq[3]
        bool have_events[3]
        bool have_model[2]
        bool have_alignment

    cppclass Cpp_File "fast5::File":
        Cpp_File() except +
        Cpp_File(string) except +
        Cpp_File(string, bool) except +

        bool is_open()
        bool is_rw()
        string file_name()
        void open(string) except +
        void open(string, bool) except +
        void create(string) except +
        void create(string, bool) except +
        void close() except +
        @staticmethod
        bool is_valid_file(string)

        string file_version() except +

        bool have_channel_id_params()
        Channel_Id_Params get_channel_id_params()
        bool have_sampling_rate()
        double get_sampling_rate()

        bool have_tracking_id_params()
        Tracking_Id_Params get_tracking_id_params() except +

        bool have_sequences_params()
        Sequences_Params get_sequences_params() except +

        vector[string] get_raw_samples_read_name_list()
        bool have_raw_samples()
        bool have_raw_samples(string)
        bool have_raw_samples_unpack(string) except +
        bool have_raw_samples_pack(string) except +
        Raw_Samples_Params get_raw_samples_params() except +
        Raw_Samples_Params get_raw_samples_params(string) except +
        vector[Raw_Int_Sample] get_raw_int_samples() except +
        vector[Raw_Int_Sample] get_raw_int_samples(string) except +
        vector[Raw_Sample] get_raw_samples() except +
        vector[Raw_Sample] get_raw_samples(string) except +

        vector[string] get_eventdetection_group_list()
        bool have_eventdetection_group()
        bool have_eventdetection_group(string)
        vector[string] get_eventdetection_read_name_list()
        vector[string] get_eventdetection_read_name_list(string)
        bool have_eventdetection_events()
        bool have_eventdetection_events(string)
        bool have_eventdetection_events(string, string)
        bool have_eventdetection_events_unpack(string, string) except +
        bool have_eventdetection_events_pack(string, string) except +
        Attr_Map get_eventdetection_params() except +
        Attr_Map get_eventdetection_params(string) except +
        EventDetection_Events_Params get_eventdetection_events_params() except +
        EventDetection_Events_Params get_eventdetection_events_params(string) except +
        EventDetection_Events_Params get_eventdetection_events_params(string, string) except +
        vector[EventDetection_Event] get_eventdetection_events() except +
        vector[EventDetection_Event] get_eventdetection_events(string) except +
        vector[EventDetection_Event] get_eventdetection_events(string, string) except +

        vector[string] get_basecall_group_list()
        bool have_basecall_group()
        bool have_basecall_group(string)
        vector[string] get_basecall_strand_group_list(unsigned)
        bool have_basecall_strand_group(unsigned)
        bool have_basecall_strand_group(unsigned, string)
        Basecall_Group_Description get_basecall_group_description(string) except +
        string get_basecall_1d_group(string)
        string get_basecall_eventdetection_group(string)
        Attr_Map get_basecall_params(string) except +
        bool have_basecall_log(string)
        string get_basecall_log(string) except +
        Attr_Map get_basecall_config(string) except +
        Attr_Map get_basecall_summary(string) except +

        bool have_basecall_fastq(unsigned)
        bool have_basecall_fastq(unsigned, string)
        bool have_basecall_fastq_unpack(unsigned, string) except +
        bool have_basecall_fastq_pack(unsigned, string) except +
        string get_basecall_fastq(unsigned) except +
        string get_basecall_fastq(unsigned, string) except +
        bool have_basecall_seq(unsigned)
        bool have_basecall_seq(unsigned, string)
        string get_basecall_seq(unsigned) except +
        string get_basecall_seq(unsigned, string) except +

        bool have_basecall_model(unsigned)
        bool have_basecall_model(unsigned, string)
        string get_basecall_model_file(unsigned) except +
        string get_basecall_model_file(unsigned, string) except +
        Basecall_Model_Params get_basecall_model_params(unsigned) except +
        Basecall_Model_Params get_basecall_model_params(unsigned, string) except +
        vector[Basecall_Model_State] get_basecall_model(unsigned) except +
        vector[Basecall_Model_State] get_basecall_model(unsigned, string) except +

        bool have_basecall_events(unsigned)
        bool have_basecall_events(unsigned, string)
        bool have_basecall_events_unpack(unsigned, string) except +
        bool have_basecall_events_pack(unsigned, string) except +
        Basecall_Events_Params get_basecall_events_params(unsigned) except +
        Basecall_Events_Params get_basecall_events_params(unsigned, string) except +
        vector[Basecall_Event] get_basecall_events(unsigned) except +
        vector[Basecall_Event] get_basecall_events(unsigned, string) except +

        bool have_basecall_alignment()
        bool have_basecall_alignment(string)
        bool have_basecall_alignment_unpack(string) except +
        bool have_basecall_alignment_pack(string) except +
        vector[Basecall_Alignment_Entry] get_basecall_alignment() except +
        vector[Basecall_Alignment_Entry] get_basecall_alignment(string) except +

cdef extern from "File_Packer.hpp" namespace "fast5":

    struct Counts "fast5::File_Packer::Counts":
        size_t rs_count
        size_t rs_bits
        size_t ed_count
        size_t ed_skip_bits
        size_t ed_len_bits
        size_t fq_count
        size_t bp_seq_count
        size_t fq_bp_bits
        size_t fq_qv_bits
        size_t ev_count
        size_t ev_rel_skip_bits
        size_t ev_skip_bits
        size_t ev_len_bits
        size_t ev_move_bits
        size_t ev_p_model_state_bits
        size_t al_count
        size_t al_template_step_bits
        size_t al_complement_step_bits
        size_t al_move_bits
        double rs_total_duration
        double rs_called_duration

    cppclass Cpp_File_Packer "fast5::File_Packer":

        Cpp_File_Packer()
        Cpp_File_Packer(int)
        Cpp_File_Packer(int, int, int, int, int)

        void set_check(bool)
        void set_force(bool)
        void set_qv_bits(unsigned)
        void set_p_model_state_bits(unsigned)

        void run(string, string) except +
        void reset_counts()
        Counts get_counts()

__version__ = cpp_version

cdef class Logger:
    @staticmethod
    def set_levels_from_options(s):
        Cpp_Logger.set_levels_from_options(s)

cdef class File:
    cdef unique_ptr[Cpp_File] thisptr

    def __init__(self, name=None, rw=None):
        if name is None:
            self.thisptr.reset(new Cpp_File())
        elif rw is None:
            self.thisptr.reset(new Cpp_File(name))
        else:
            self.thisptr.reset(new Cpp_File(name, rw))

    def is_open(self):
        return deref(self.thisptr).is_open()
    def is_rw(self):
        return deref(self.thisptr).is_rw()
    def file_name(self):
        return deref(self.thisptr).file_name()
    def open(self, file_name, rw=None):
        if rw is None:
            return deref(self.thisptr).open(file_name)
        else:
            return deref(self.thisptr).open(file_name, rw)
    def create(self, file_name, trunc=None):
        if trunc is None:
            return deref(self.thisptr).open(file_name)
        else:
            return deref(self.thisptr).open(file_name, trunc)
    def close(self):
        return deref(self.thisptr).close()
    @staticmethod
    def is_valid_file(s):
        return Cpp_File.is_valid_file(s)

    def file_version(self):
        return deref(self.thisptr).file_version()

    def have_channel_id_params(self):
        return deref(self.thisptr).have_channel_id_params()
    def get_channel_id_params(self):
        return deref(self.thisptr).get_channel_id_params()

    def have_tracking_id_params(self):
        return deref(self.thisptr).have_tracking_id_params()
    def get_tracking_id_params(self):
        return deref(self.thisptr).get_tracking_id_params()

    def have_sequences_params(self):
        return deref(self.thisptr).have_sequences_params()
    def get_sequences_params(self):
        return deref(self.thisptr).get_sequences_params()

    def get_raw_samples_read_name_list(self):
        return deref(self.thisptr).get_raw_samples_read_name_list()
    def have_raw_samples(self, rn=None):
        if rn is None:
            return deref(self.thisptr).have_raw_samples()
        else:
            return deref(self.thisptr).have_raw_samples(rn)
    def have_raw_samples_unpack(self, rn):
        return deref(self.thisptr).have_raw_samples_unpack(rn)
    def have_raw_samples_pack(self, rn):
        return deref(self.thisptr).have_raw_samples_pack(rn)
    def get_raw_samples_params(self, rn=None):
        if rn is None:
            return deref(self.thisptr).get_raw_samples_params()
        else:
            return deref(self.thisptr).get_raw_samples_params(rn)
    def get_raw_int_samples(self, rn=None):
        if rn is None:
            return deref(self.thisptr).get_raw_int_samples()
        else:
            return deref(self.thisptr).get_raw_int_samples(rn)
    def get_raw_samples(self, rn=None):
        if rn is None:
            return deref(self.thisptr).get_raw_samples()
        else:
            return deref(self.thisptr).get_raw_samples(rn)

    def get_eventdetection_group_list(self):
        return deref(self.thisptr).get_eventdetection_group_list()
    def have_eventdetection_group(self, gr=None):
        if gr is None:
            return deref(self.thisptr).have_eventdetection_group()
        else:
            return deref(self.thisptr).have_eventdetection_group(gr)
    def get_eventdetection_params(self, gr=None):
        if gr is None:
            return deref(self.thisptr).get_eventdetection_params()
        else:
            return deref(self.thisptr).get_eventdetection_params(gr)
    def get_eventdetection_read_name_list(self, gr=None):
        if gr is None:
            return deref(self.thisptr).get_eventdetection_read_name_list()
        else:
            return deref(self.thisptr).get_eventdetection_read_name_list(gr)
    def have_eventdetection_events(self, gr=None, rn=None):
        if gr is None:
            return deref(self.thisptr).have_eventdetection_events()
        elif rn is None:
            return deref(self.thisptr).have_eventdetection_events(gr)
        else:
            return deref(self.thisptr).have_eventdetection_events(gr, rn)
    def have_eventdetection_events_unpack(self, gr, rn):
        return deref(self.thisptr).have_eventdetection_events_unpack(gr, rn)
    def have_eventdetection_events_pack(self, gr, rn):
        return deref(self.thisptr).have_eventdetection_events_pack(gr, rn)
    def get_eventdetection_events_params(self, gr=None, rn=None):
        if gr is None:
            return deref(self.thisptr).get_eventdetection_events_params()
        elif rn is None:
            return deref(self.thisptr).get_eventdetection_events_params(gr)
        else:
            return deref(self.thisptr).get_eventdetection_events_params(gr, rn)
    def get_eventdetection_events(self, gr=None, rn=None):
        if gr is None:
            return deref(self.thisptr).get_eventdetection_events()
        elif rn is None:
            return deref(self.thisptr).get_eventdetection_events(gr)
        else:
            return deref(self.thisptr).get_eventdetection_events(gr, rn)

    def get_basecall_group_list(self):
        return deref(self.thisptr).get_basecall_group_list()
    def get_basecall_strand_group_list(self, st):
        return deref(self.thisptr).get_basecall_strand_group_list(st)
    def have_basecall_group(self, gr=None):
        if gr is None:
            return deref(self.thisptr).have_basecall_group()
        else:
            return deref(self.thisptr).have_basecall_group(gr)
    def have_basecall_strand_group(self, st, gr=None):
        if gr is None:
            return deref(self.thisptr).have_basecall_strand_group(st)
        else:
            return deref(self.thisptr).have_basecall_strand_group(st, gr)
    def get_basecall_group_description(self, gr):
        return deref(self.thisptr).get_basecall_group_description(gr)
    def get_basecall_1d_group(self, gr):
        return deref(self.thisptr).get_basecall_1d_group(gr)
    def get_basecall_eventdetection_group(self, gr):
        return deref(self.thisptr).get_basecall_eventdetection_group(gr)
    def get_basecall_params(self, gr):
        return deref(self.thisptr).get_basecall_params(gr)
    def get_basecall_log(self, gr):
        return deref(self.thisptr).get_basecall_log(gr)
    def get_basecall_config(self, gr):
        return deref(self.thisptr).get_basecall_config(gr)
    def get_basecall_summary(self, gr):
        return deref(self.thisptr).get_basecall_summary(gr)

    def have_basecall_fastq(self, st, gr=None):
        if gr is None:
            return deref(self.thisptr).have_basecall_fastq(st)
        else:
            return deref(self.thisptr).have_basecall_fastq(st, gr)
    def have_basecall_fastq_unpack(self, st, gr):
        return deref(self.thisptr).have_basecall_fastq_unpack(st, gr)
    def have_basecall_fastq_pack(self, st, gr):
        return deref(self.thisptr).have_basecall_fastq_pack(st, gr)
    def get_basecall_fastq(self, st, gr=None):
        if gr is None:
            return deref(self.thisptr).get_basecall_fastq(st)
        else:
            return deref(self.thisptr).get_basecall_fastq(st, gr)
    def have_basecall_seq(self, st, gr=None):
        if gr is None:
            return deref(self.thisptr).have_basecall_seq(st)
        else:
            return deref(self.thisptr).have_basecall_seq(st, gr)
    def get_basecall_seq(self, st, gr=None):
        if gr is None:
            return deref(self.thisptr).get_basecall_seq(st)
        else:
            return deref(self.thisptr).get_basecall_seq(st, gr)

    def have_basecall_model(self, st, gr=None):
        if gr is None:
            return deref(self.thisptr).have_basecall_model(st)
        else:
            return deref(self.thisptr).have_basecall_model(st, gr)
    def get_basecall_model_file(self, st, gr=None):
        if gr is None:
            return deref(self.thisptr).get_basecall_model_file(st)
        else:
            return deref(self.thisptr).get_basecall_model_file(st, gr)
    def get_basecall_model_params(self, st, gr=None):
        if gr is None:
            return deref(self.thisptr).get_basecall_model_params(st)
        else:
            return deref(self.thisptr).get_basecall_model_params(st, gr)
    def get_basecall_model(self, st, gr=None):
        if gr is None:
            return deref(self.thisptr).get_basecall_model(st)
        else:
            return deref(self.thisptr).get_basecall_model(st, gr)

    def have_basecall_events(self, st, gr=None):
        if gr is None:
            return deref(self.thisptr).have_basecall_events(st)
        else:
            return deref(self.thisptr).have_basecall_events(st, gr)
    def have_basecall_events_unpack(self, st, gr):
        return deref(self.thisptr).have_basecall_events_unpack(st, gr)
    def have_basecall_events_pack(self, st, gr):
        return deref(self.thisptr).have_basecall_events_pack(st, gr)
    def get_basecall_events_params(self, st, gr=None):
        if gr is None:
            return deref(self.thisptr).get_basecall_events_params(st)
        else:
            return deref(self.thisptr).get_basecall_events_params(st, gr)
    def get_basecall_events(self, st, gr=None):
        if gr is None:
            return deref(self.thisptr).get_basecall_events(st)
        else:
            return deref(self.thisptr).get_basecall_events(st, gr)

    def have_basecall_alignment(self, gr=None):
        if gr is None:
            return deref(self.thisptr).have_basecall_alignment()
        else:
            return deref(self.thisptr).have_basecall_alignment(gr)
    def have_basecall_alignment_unpack(self, gr):
        return deref(self.thisptr).have_basecall_alignment_unpack(gr)
    def have_basecall_alignment_pack(self, gr):
        return deref(self.thisptr).have_basecall_alignment_pack(gr)
    def get_basecall_alignment(self, gr=None):
        if gr is None:
            return deref(self.thisptr).get_basecall_alignment()
        else:
            return deref(self.thisptr).get_basecall_alignment(gr)

cdef class File_Packer:
    cdef unique_ptr[Cpp_File_Packer] thisptr

    def __init__(self, a1=None, a2=None, a3=None, a4=None, a5=None):
        if a1 is None:
            self.thisptr.reset(new Cpp_File_Packer())
        elif a2 is None:
            self.thisptr.reset(new Cpp_File_Packer(a1))
        else:
            self.thisptr.reset(new Cpp_File_Packer(a1, a2, a3, a4, a5))

    def set_check(self, _check):
        deref(self.thisptr).set_check(_check)
    def set_force(self, _force):
        deref(self.thisptr).set_force(_force)
    def set_qv_bits(self, _qv_bits):
        deref(self.thisptr).set_qv_bits(_qv_bits)
    def set_p_model_state_bits(self, _p_model_state_bits):
        deref(self.thisptr).set_p_model_state_bits(_p_model_state_bits)

    def run(self, ifn, ofn):
        deref(self.thisptr).run(ifn, ofn)
    def reset_counts(self):
        deref(self.thisptr).reset_counts()
    def get_counts(self):
        return deref(self.thisptr).get_counts()
