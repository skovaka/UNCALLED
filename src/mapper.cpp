/* MIT License
 *
 * Copyright (c) 2018 Sam Kovaka <skovaka@gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "pdqsort.h"
#include "mapper.hpp"

MapperParams::MapperParams(const std::string &bwa_prefix,
                           const std::string &model_fname,
                           u32 seed_len, 
                           u32 min_aln_len,
                           u32 min_rep_len, 
                           u32 max_rep_copy, 
                           u32 max_consec_stay,
                           u32 max_paths, 
                           u32 max_events_proc,
                           u32 evt_buffer_len,
                           u32 evt_winlen1,
                           u32 evt_winlen2,
                           u16 evt_batch_size,
                           float evt_timeout,
                           float evt_thresh1,
                           float evt_thresh2,
                           float evt_peak_height,
                           float evt_min_mean,
                           float evt_max_mean,
                           float max_stay_frac,
                           float min_seed_prob, 
                           float min_mean_conf,
                           float min_top_conf)
        : fmi_(BwaFMI(bwa_prefix)),
          model_(KmerModel(model_fname, true)),
          event_params_({evt_winlen1,
                         evt_winlen2,
                         evt_thresh1,
                         evt_thresh2,
                         evt_peak_height,
                         evt_min_mean,
                         evt_max_mean}),
          seed_len_(seed_len),
          min_rep_len_(min_rep_len),
          max_rep_copy_(max_rep_copy),
          max_paths_(max_paths),
          max_consec_stay_(max_consec_stay),
          min_aln_len_(min_aln_len),
          max_events_proc_(max_events_proc),
          evt_buffer_len_(evt_buffer_len),
          evt_batch_size_(evt_batch_size),
          evt_timeout_(evt_timeout),
          max_stay_frac_(max_stay_frac),
          min_seed_prob_(min_seed_prob),
          min_mean_conf_(min_mean_conf),
          min_top_conf_(min_top_conf) {
    
    //TODO: exception handling
    std::ifstream infile(bwa_prefix + INDEX_SUFF);
    float prob, frac;
    u64 fmlen = 0;
    infile >> prob >> frac;
    evpr_threshes_.push_back(prob);
    while (fmlen != 1) {
        infile >> fmlen >> prob >> frac;
        evpr_lengths_.push_back(fmlen);
        evpr_threshes_.push_back(prob);
    }

    kmer_fmranges_ = std::vector<Range>(model_.kmer_count());
    for (u16 k = 0; k < model_.kmer_count(); k++) {
        Range r = fmi_.get_full_range(model_.get_last_base(k));
        for (u8 i = model_.kmer_len()-2; i < model_.kmer_len(); i--) {
            r = fmi_.get_neighbor(r, model_.get_base(k, i));
        }
        kmer_fmranges_[k] = r;
    }
}

float MapperParams::get_prob_thresh(u64 fm_length) const {
    auto pr = evpr_threshes_.begin();
    for (auto len = evpr_lengths_.begin(); len != evpr_lengths_.end(); len++) {
        if (fm_length > *len) {
            break;
        }
        pr++;
    }
    return *pr;
}

float MapperParams::get_source_prob() const {
    return evpr_threshes_.front();
}

u16 MapperParams::get_max_events(u16 event_i) const {
    if (event_i + evt_batch_size_ > max_events_proc_) 
        return max_events_proc_ - event_i;
    return evt_batch_size_;
}

ReadLoc::ReadLoc(const std::string &rd_name, u16 channel, u32 number) 
    : rd_name_(rd_name),
      rd_channel_(channel),
      rd_number_(number),
      rd_st_(0),
      rd_en_(0),
      rd_len_(0),
      match_count_(0),
      time_(-1) {
    
    #ifdef DEBUG_TIME
    sigproc_time_ = prob_time_ = thresh_time_ = stay_time_ = 
    neighbor_time_ = fmrs_time_ = fmsa_time_ = 
    sort_time_ = loop2_time_ = source_time_ = tracker_time_ = 0;
    #endif
}

ReadLoc::ReadLoc() {
    match_count_ = 0;

    #ifdef DEBUG_TIME
    sigproc_time_ = prob_time_ = thresh_time_ = stay_time_ = 
    neighbor_time_ = fmrs_time_ = fmsa_time_ = 
    sort_time_ = loop2_time_ = source_time_ = tracker_time_ = 0;
    #endif
}

bool ReadLoc::set_ref_loc(const MapperParams &params, const SeedGroup &seeds) {
    u8 k_shift = (params.model_.kmer_len() - 1);

    rd_st_ = (u32) (params.max_stay_frac_ * seeds.evt_st_);
    rd_en_ = (u32) (params.max_stay_frac_ * (seeds.evt_en_ + params.seed_len_)) + k_shift;

    match_count_ = seeds.total_len_ + k_shift;
    fwd_ = seeds.ref_st_ > params.fmi_.size() / 2;

    u64 sa_st;
    if (fwd_) sa_st = params.fmi_.size() - (seeds.ref_en_.end_ + k_shift);
    else      sa_st = seeds.ref_st_;

    rf_len_ = params.fmi_.translate_loc(sa_st, rf_name_, rf_st_);
    rf_en_ = rf_st_ + (seeds.ref_en_.end_ - seeds.ref_st_) + k_shift;

    return rf_len_ > 0;
}

void ReadLoc::set_time(float time) {
    time_ = time;
}

void ReadLoc::set_read_len(const MapperParams &params, u32 len) {
    rd_len_ = (u32) (params.max_stay_frac_ * len);
    rd_en_ = rd_len_ - 1;
}

bool ReadLoc::is_valid() const {
    return match_count_ > 0;
}

u16 ReadLoc::get_channel() const {
    return rd_channel_;
}

u32 ReadLoc::get_number() const {
    return rd_number_;
}

std::string ReadLoc::str() const {
    std::stringstream ss;
    ss << rd_name_ << "\t"
       << rd_len_ << "\t"
       << rd_st_ << "\t"
       << rd_en_ << "\t";
    if (is_valid()) {
        ss << (fwd_ ? '+' : '-') << "\t"
           << rf_name_ << "\t"
           << rf_len_ << "\t"
           << rf_st_ << "\t"
           << rf_en_ << "\t"
           << match_count_ << "\t"
           << (rf_en_ - rf_st_ + 1) << "\t"
           << 255;
    } else {
        ss << "*" << "\t"
           << "*" << "\t"
           << "*" << "\t"
           << "*" << "\t"
           << "*" << "\t"
           << "*" << "\t"
           << "*" << "\t"
           << "255";
    }

    #ifndef DEBUG_TIME
    if (time_ > 0) {
        ss << "\t" << PAF_TIME_TAG << time_;
    }
    #else
        ss << "\t" PAF_SIGPROC_TAG << sigproc_time_  
           << "\t" PAF_PROB_TAG    << prob_time_
           << "\t" PAF_THRESH_TAG    << thresh_time_
           << "\t" PAF_STAY_TAG    << stay_time_
           << "\t" PAF_NEIGHBOR_TAG    << neighbor_time_
           << "\t" PAF_FMRS_TAG    << fmrs_time_
           << "\t" PAF_FMSA_TAG    << fmsa_time_
           << "\t" PAF_SORT_TAG    << sort_time_
           << "\t" PAF_LOOP2_TAG   << loop2_time_
           << "\t" PAF_SOURCE_TAG  << source_time_
           << "\t" PAF_TRACKER_TAG << tracker_time_;
    #endif

    return ss.str();
}

std::ostream &operator<< (std::ostream &out, const ReadLoc &l) {
    out << l.str();
    return out;
}

u8 Mapper::PathBuffer::MAX_PATH_LEN = 0, 
   Mapper::PathBuffer::TYPE_MASK = 0;

u32 Mapper::PathBuffer::TYPE_ADDS[EventType::NUM_TYPES];

Mapper::PathBuffer::PathBuffer()
    : length_(0),
      prob_sums_(new float[MAX_PATH_LEN+1]) {
}

Mapper::PathBuffer::PathBuffer(const PathBuffer &p) {
    std::memcpy(this, &p, sizeof(PathBuffer));
}

void Mapper::PathBuffer::free_buffers() {
    delete[] prob_sums_;
}

void Mapper::PathBuffer::make_source(Range &range, u16 kmer, float prob) {
    length_ = 1;
    consec_stays_ = 0;
    event_types_ = 0;
    seed_prob_ = prob;
    fm_range_ = range;
    kmer_ = kmer;
    sa_checked_ = false;

    path_type_counts_[EventType::MATCH] = 1;

    //TODO: no loops!
    for (u8 t = 1; t < EventType::NUM_TYPES; t++) {
        path_type_counts_[t] = 0;
    }

    //TODO: don't write this here to speed up source loop
    prob_sums_[0] = 0;
    prob_sums_[1] = prob;
}


void Mapper::PathBuffer::make_child(PathBuffer &p, 
                                    Range &range,
                                    u16 kmer, 
                                    float prob, 
                                    EventType type) {

    length_ = p.length_ + (p.length_ <= MAX_PATH_LEN);
    fm_range_ = range;
    kmer_ = kmer;
    sa_checked_ = p.sa_checked_;
    event_types_ = TYPE_ADDS[type] | (p.event_types_ >> TYPE_BITS);
    consec_stays_ = (p.consec_stays_ + (type == EventType::STAY)) * (type == EventType::STAY);

    std::memcpy(path_type_counts_, p.path_type_counts_, EventType::NUM_TYPES * sizeof(u8));
    path_type_counts_[type]++;


    if (length_ > MAX_PATH_LEN) {
        std::memcpy(prob_sums_, &(p.prob_sums_[1]), MAX_PATH_LEN * sizeof(float));
        prob_sums_[MAX_PATH_LEN] = prob_sums_[MAX_PATH_LEN-1] + prob;
        seed_prob_ = (prob_sums_[MAX_PATH_LEN] - prob_sums_[0]) / MAX_PATH_LEN;
        path_type_counts_[p.type_tail()]--;
    } else {
        std::memcpy(prob_sums_, p.prob_sums_, length_ * sizeof(float));
        prob_sums_[length_] = prob_sums_[length_-1] + prob;
        seed_prob_ = prob_sums_[length_] / length_;
    }


}

void Mapper::PathBuffer::invalidate() {
    length_ = 0;
}

bool Mapper::PathBuffer::is_valid() const {
    return length_ > 0;
}

u8 Mapper::PathBuffer::match_len() const {
    return path_type_counts_[EventType::MATCH];
}

u8 Mapper::PathBuffer::type_head() const {
    return (event_types_ >> (TYPE_BITS*(MAX_PATH_LEN-2))) & TYPE_MASK;
}

u8 Mapper::PathBuffer::type_tail() const {
    return event_types_ & TYPE_MASK;
}

bool Mapper::PathBuffer::is_seed_valid(const MapperParams &p,
                                        bool path_ended) const{
    return (fm_range_.length() == 1 || 
                (path_ended &&
                 fm_range_.length() <= p.max_rep_copy_ &&
                 match_len() >= p.min_rep_len_)) &&

           length_ >= p.seed_len_ &&
           (path_ended || type_head() == EventType::MATCH) &&
           (path_ended || path_type_counts_[EventType::STAY] <= p.max_stay_frac_ * p.seed_len_) &&
          seed_prob_ >= p.min_seed_prob_;
}


bool operator< (const Mapper::PathBuffer &p1, 
                const Mapper::PathBuffer &p2) {
    return p1.fm_range_ < p2.fm_range_ ||
           (p1.fm_range_ == p2.fm_range_ && 
            p1.seed_prob_ < p2.seed_prob_);
}

Mapper::Mapper(const MapperParams &ap, u16 channel)
    : 
      params_(ap),
      model_(ap.model_),
      fmi_(ap.fmi_),
      event_detector_(ap.event_params_),
      norm_(ap.model_, ap.evt_buffer_len_),
      seed_tracker_(ap.fmi_.size(),
                    ap.min_mean_conf_,
                    ap.min_top_conf_,
                    ap.min_aln_len_,
                    ap.seed_len_),
      channel_(channel),
      read_num_(0),
      chunk_processed_(true),
      state_(State::INACTIVE)

     {


    PathBuffer::MAX_PATH_LEN = params_.seed_len_;

    for (u64 t = 0; t < EventType::NUM_TYPES; t++) {
        PathBuffer::TYPE_ADDS[t] = t << ((PathBuffer::MAX_PATH_LEN-2)*TYPE_BITS);
    }
    PathBuffer::TYPE_MASK = (u8) ((1 << TYPE_BITS) - 1);

    kmer_probs_ = std::vector<float>(model_.kmer_count());
    prev_paths_ = std::vector<PathBuffer>(params_.max_paths_);
    next_paths_ = std::vector<PathBuffer>(params_.max_paths_);
    sources_added_ = std::vector<bool>(model_.kmer_count(), false);

    prev_size_ = 0;
    event_i_ = 0;
    seed_tracker_.reset();
}

Mapper::Mapper(const Mapper &m) : Mapper(m.params_, m.channel_) {}

Mapper::~Mapper() {
    for (u32 i = 0; i < next_paths_.size(); i++) {
        next_paths_[i].free_buffers();
        prev_paths_[i].free_buffers();
    }
}

void Mapper::skip_events(u32 n) {
    event_i_ += n;
    prev_size_ = 0;
}

void Mapper::new_read(const std::string &id, u32 number) {
    read_loc_ = ReadLoc(id, channel_, number);
    read_num_ = number;
    prev_size_ = 0;
    event_i_ = 0;
    chunk_processed_ = true;
    reset_ = false;
    last_chunk_ = false;
    state_ = State::MAPPING;
    seed_tracker_.reset();
    event_detector_.reset();
    norm_.skip_unread();
    chunk_buffer_.clear();
    timer_.reset();
}

std::string Mapper::map_fast5(const std::string &fast5_name) {
    if (!fast5::File::is_valid_file(fast5_name)) {
        std::cerr << "Error: '" << fast5_name << "' is not a valid file \n";
    }

    ReadLoc aln;

    try {
        fast5::File file;
        file.open(fast5_name);
        
        if (file.is_open()) {  
            auto fast5_info = file.get_raw_samples_params();
            auto raw_samples = file.get_raw_samples();
            new_read(fast5_info.read_id, fast5_info.read_number);
            aln = add_samples(raw_samples);

        } else {
            std::cerr << "Error: unable to open '" << fast5_name << "'\n";
        }

        
    } catch (hdf5_tools::Exception& e) {
        std::cerr << "Error: hdf5 exception '" << e.what() << "'\n";
    }

    return aln.str();
}

bool Mapper::add_sample(float s) {
    ;
    if (!event_detector_.add_sample(s)) return false;
    norm_.add_event(event_detector_.get_mean());
    float m = norm_.pop_event();

    #ifdef DEBUG_TIME
    read_loc_.sigproc_time_ += timer_.lap();
    #endif

    if (event_i_ >= params_.max_events_proc_ || add_event(m)) {
        read_loc_.set_time(timer_.get());
        read_loc_.set_read_len(params_, event_i_);
        return true;
    }

    return false;
}

u32 Mapper::prev_unfinished(u32 next_number) const {
    return state_ == State::MAPPING && read_num_ != next_number;
}

bool Mapper::finished() const {
    return state_ == State::SUCCESS || state_ == State::FAILURE;
}

ReadLoc Mapper::pop_loc() {
    state_ = State::INACTIVE;
    reset_ = false;
    return read_loc_;
}

ReadLoc Mapper::get_loc() const {
    return read_loc_;
}

ReadLoc Mapper::add_samples(const std::vector<float> &samples) {

    if (params_.evt_buffer_len_ == 0) {
        #ifdef DEBUG_TIME
        timer_.reset();
        #endif

        std::vector<Event> events = event_detector_.add_samples(samples);
        std::vector<Event> old(events);
        model_.normalize(events);

        #ifdef DEBUG_TIME
        read_loc_.sigproc_time_ += timer_.lap();
        #endif

        read_loc_.set_read_len(params_, events.size());

        for (u32 e = 0; e < events.size(); e++) {
            if (add_event(events[e].mean)) break;
        }
    } else {

        read_loc_.set_read_len(params_, samples.size() / 5); //TODO: this better
        
        u32 i = 0;

        #ifdef DEBUG_TIME
        timer_.reset();
        #endif

        float m;
        for (auto s : samples) {
            if (!event_detector_.add_sample(s)) continue;
            norm_.add_event(event_detector_.get_mean());
            m = norm_.pop_event();

            #ifdef DEBUG_TIME
            read_loc_.sigproc_time_ += timer_.lap();
            #endif

            if (add_event(m)) break;
            i++;
        }
    }

    read_loc_.set_time(timer_.get());

    return read_loc_;
}

void Mapper::request_reset() {
    reset_ = true;
}

void Mapper::end_reset() {
    reset_ = false;
}

bool Mapper::is_resetting() {
    return reset_;
}

bool Mapper::is_chunk_processed() const {
    return chunk_processed_;
}

Mapper::State Mapper::get_state() const {
    return state_;
}

bool Mapper::swap_chunk(Chunk &chunk) {
    //std::cout << "# tryna swap " << is_chunk_processed() << " " << reset_ << " " << state_ << "\n";
    if (!is_chunk_processed() || reset_) return false;

    //New read hasn't been set
    if (chunk.number != read_num_) {
        state_ = State::FAILURE;
    }

    chunk_buffer_.swap(chunk.raw_data);
    chunk_processed_ = false;
    return true;
}

u16 Mapper::process_chunk() {
    if (chunk_processed_ || reset_) return 0; 
    
    #ifdef DEBUG_TIME
    read_loc_.sigproc_time_ += timer_.lap();
    #endif

    float mean;

    //std::cout << "# processing " << channel_ << "\n";

    u16 nevents = 0;
    for (u32 i = 0; i < chunk_buffer_.size(); i++) {
        if (event_detector_.add_sample(chunk_buffer_[i])) {
            mean = event_detector_.get_mean();
            if (!norm_.add_event(mean)) {

                //std::cout << "# skipping unread " << channel_ << "\n";
                u32 nskip = norm_.skip_unread(nevents);
                skip_events(nskip);
                //TODO: report event skip in some way
                //std::cout << "# norm skipped " << nskip << "\n";
                if (!norm_.add_event(mean)) {
                    std::cout << "# error: chunk events cannot fit in normilzation buffer\n";
                    return nevents;
                }
            }
            nevents++;
        }
    }

    //std::cout << "# processed " << channel_ << "\n";

    chunk_buffer_.clear();
    chunk_processed_ = true;
    return nevents;
}

bool Mapper::end_read(u32 number) {
    //set last chunk if you want to keep trying after read has ended
    //return last_chunk_ = (read_loc_.get_number() == number);
    return reset_ = (read_loc_.get_number() == number);
}

bool Mapper::map_chunk() {
    if (reset_ || (last_chunk_ && norm_.empty())) return true;
    u16 nevents = params_.get_max_events(event_i_);
    float tlimit = params_.evt_timeout_ * nevents;

    Timer t;
    for (u16 i = 0; i < nevents && !norm_.empty(); i++) {
        if (add_event(norm_.pop_event())) return true;
        if (t.get() > tlimit) {
            //std::cout << "# timeout " << channel_ << " " << i << "\n";
            return false; //TODO: penalize this read
        }
    }

    return false;
}

bool Mapper::add_event(float event) {

    if (reset_ || event_i_ >= params_.max_events_proc_) {
        reset_ = false;
        state_ = State::FAILURE;
        return true;
    }

    Range prev_range;
    u16 prev_kmer;
    float evpr_thresh;
    bool child_found;


    auto next_path = next_paths_.begin();

    for (u16 kmer = 0; kmer < model_.kmer_count(); kmer++) {
        kmer_probs_[kmer] = model_.event_match_prob(event, kmer);
    }
    #ifdef DEBUG_TIME
    read_loc_.prob_time_ += timer_.lap();
    #endif
    
    //Find neighbors of previous nodes
    for (u32 pi = 0; pi < prev_size_; pi++) {
        if (!prev_paths_[pi].is_valid()) {
            continue;
        }

        child_found = false;

        PathBuffer &prev_path = prev_paths_[pi];
        Range &prev_range = prev_path.fm_range_;
        prev_kmer = prev_path.kmer_;

        evpr_thresh = params_.get_prob_thresh(prev_range.length());
        #ifdef DEBUG_TIME
        read_loc_.thresh_time_ += timer_.lap();
        #endif

        if (prev_path.consec_stays_ < params_.max_consec_stay_ && 
            kmer_probs_[prev_kmer] >= evpr_thresh) {

            next_path->make_child(prev_path, 
                                  prev_range,
                                  prev_kmer, 
                                  kmer_probs_[prev_kmer], 
                                  EventType::STAY);

            child_found = true;

            if (++next_path == next_paths_.end()) {
                break;
            }
        }


        #ifdef DEBUG_TIME
        read_loc_.stay_time_ += timer_.lap();
        #endif

        //Add all the neighbors
        for (u8 b = 0; b < ALPH_SIZE; b++) {
            u16 next_kmer = model_.get_neighbor(prev_kmer, b);

            if (kmer_probs_[next_kmer] < evpr_thresh) {
                continue;
            }

            #ifdef DEBUG_TIME
            read_loc_.neighbor_time_ += timer_.lap();
            #endif

            Range next_range = fmi_.get_neighbor(prev_range, b);

            #ifdef DEBUG_TIME
            read_loc_.fmrs_time_ += timer_.lap();
            #endif

            if (!next_range.is_valid()) {
                continue;
            }

            next_path->make_child(prev_path, 
                                  next_range,
                                  next_kmer, 
                                  kmer_probs_[next_kmer], 
                                  EventType::MATCH);

            child_found = true;

            if (++next_path == next_paths_.end()) {
                break;
            }
        }

        #ifdef DEBUG_TIME
        read_loc_.neighbor_time_ += timer_.lap();
        #endif

        if (!child_found && !prev_path.sa_checked_) {

            update_seeds(prev_path, true);

        }

        if (next_path == next_paths_.end()) {
            break;
        }
    }

    if (next_path != next_paths_.begin()) {

        u32 next_size = next_path - next_paths_.begin();

        pdqsort(next_paths_.begin(), next_path);
        //std::sort(next_paths_.begin(), next_path);

        #ifdef DEBUG_TIME
        read_loc_.sort_time_ += timer_.lap();
        #endif

        u16 source_kmer;
        prev_kmer = model_.kmer_count(); 

        Range unchecked_range, source_range;

        for (u32 i = 0; i < next_size; i++) {
            source_kmer = next_paths_[i].kmer_;

            //Add source for beginning of kmer range
            if (source_kmer != prev_kmer &&
                next_path != next_paths_.end() &&
                kmer_probs_[source_kmer] >= params_.get_source_prob()) {

                sources_added_[source_kmer] = true;

                source_range = Range(params_.kmer_fmranges_[source_kmer].start_,
                                     next_paths_[i].fm_range_.start_ - 1);

                if (source_range.is_valid()) {
                    next_path->make_source(source_range,
                                           source_kmer,
                                           kmer_probs_[source_kmer]);
                    next_path++;
                }                                    

                unchecked_range = Range(next_paths_[i].fm_range_.end_ + 1,
                                        params_.kmer_fmranges_[source_kmer].end_);
            }

            prev_kmer = source_kmer;

            //Range next_range = next_paths_[i].fm_range_;

            //Remove paths with duplicate ranges
            //Best path will be listed last
            if (i < next_size - 1 && next_paths_[i].fm_range_ == next_paths_[i+1].fm_range_) {
                next_paths_[i].invalidate();
                continue;
            }

            //Start source after current path
            //TODO: check if theres space for a source here, instead of after extra work?
            if (next_path != next_paths_.end() &&
                kmer_probs_[source_kmer] >= params_.get_source_prob()) {
                
                source_range = unchecked_range;
                
                //Between this and next path ranges
                if (i < next_size - 1 && source_kmer == next_paths_[i+1].kmer_) {

                    source_range.end_ = next_paths_[i+1].fm_range_.start_ - 1;

                    if (unchecked_range.start_ <= next_paths_[i+1].fm_range_.end_) {
                        unchecked_range.start_ = next_paths_[i+1].fm_range_.end_ + 1;
                    }
                }

                //Add it if it's a real range
                if (source_range.is_valid()) {
                    next_path->make_source(source_range,
                                           source_kmer,
                                           kmer_probs_[source_kmer]);
                    next_path++;
                }
            }

            #ifdef DEBUG_TIME
            read_loc_.loop2_time_ += timer_.lap();
            #endif

            update_seeds(next_paths_[i], false);

        }
    }

    #ifdef DEBUG_TIME
    read_loc_.loop2_time_ += timer_.lap();
    #endif
    
    for (u16 kmer = 0; 
         kmer < model_.kmer_count() && 
            next_path != next_paths_.end(); 
         kmer++) {

        Range next_range = params_.kmer_fmranges_[kmer];

        if (!sources_added_[kmer] && 
            kmer_probs_[kmer] >= params_.get_source_prob() &&
            next_path != next_paths_.end() &&
            next_range.is_valid()) {

            //TODO: don't write to prob buffer here to speed up source loop
            next_path->make_source(next_range, kmer, kmer_probs_[kmer]);
            next_path++;

        } else {
            sources_added_[kmer] = false;
        }
    }

    #ifdef DEBUG_TIME
    read_loc_.source_time_ += timer_.lap();
    #endif

    prev_size_ = next_path - next_paths_.begin();
    prev_paths_.swap(next_paths_);

    //Update event index
    event_i_++;

    SeedGroup sg = seed_tracker_.get_final();

    if (sg.is_valid()) {
        state_ = State::SUCCESS;
        read_loc_.set_ref_loc(params_, sg);
        read_num_ = 0;

        #ifdef DEBUG_TIME
        read_loc_.tracker_time_ += timer_.lap();
        #endif

        return true;
    }

    #ifdef DEBUG_TIME
    read_loc_.tracker_time_ += timer_.lap();
    #endif

    return false;
}

void Mapper::update_seeds(PathBuffer &p, bool path_ended) {

    if (p.is_seed_valid(params_, path_ended)) {

        #ifdef DEBUG_TIME
        read_loc_.tracker_time_ += timer_.lap();
        #endif

        p.sa_checked_ = true;

        for (u64 s = p.fm_range_.start_; s <= p.fm_range_.end_; s++) {

            //Reverse the reference coords so they both go L->R
            u64 ref_en = fmi_.size() - fmi_.sa(s) + 1;

            #ifdef DEBUG_TIME
            read_loc_.fmsa_time_ += timer_.lap();
            #endif

            seed_tracker_.add_seed(ref_en, p.match_len(), event_i_ - path_ended);

            #ifdef DEBUG_TIME
            read_loc_.tracker_time_ += timer_.lap();
            #endif

            #ifdef DEBUG_SEEDS
            seed.print(seeds_out);
            #endif
        }
    }

}


