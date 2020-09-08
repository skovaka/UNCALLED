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

#include <pdqsort.h>
#include "mapper.hpp"
#include "model_r94.inl"

Mapper::Params Mapper::PRMS = {
    seed_len        : 22,
    min_rep_len     : 0,
    max_rep_copy    : 50,
    max_paths       : 10000,
    max_consec_stay : 8,
    max_events      : 30000,
    max_stay_frac   : 0.5,
    min_seed_prob   : -3.75,
    evt_buffer_len  : 6000,
    evt_batch_size  : 5,
    evt_timeout     : 10.0,
    chunk_timeout  : 4000.0 ,
    bwa_prefix      : "",
    idx_preset      : "default",
    seed_prms       : SeedTracker::PRMS_DEF,
    event_prms      : EventDetector::PRMS_DEF
};

BwaIndex<KLEN> Mapper::fmi;
std::vector<float> Mapper::prob_threshes_;

PoreModel<KLEN> Mapper::model = pmodel_r94_complement;

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
    path_type_counts_[EventType::STAY] = 0;
    total_match_len_ = 1;

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
    total_match_len_ = p.total_match_len_ + (type==EventType::MATCH);

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

bool Mapper::PathBuffer::is_seed_valid(bool path_ended) const{
    return (fm_range_.length() == 1 || 
                (path_ended &&
                 fm_range_.length() <= PRMS.max_rep_copy &&
                 match_len() >= PRMS.min_rep_len)) &&

           length_ >= PRMS.seed_len &&
           (path_ended || type_head() == EventType::MATCH) &&
           (path_ended || path_type_counts_[EventType::STAY] <= PRMS.max_stay_frac * PRMS.seed_len) &&
          seed_prob_ >= PRMS.min_seed_prob;
}


bool operator< (const Mapper::PathBuffer &p1, 
                const Mapper::PathBuffer &p2) {
    return p1.fm_range_ < p2.fm_range_ ||
           (p1.fm_range_ == p2.fm_range_ && 
            p1.seed_prob_ < p2.seed_prob_);
}

Mapper::Mapper() :
    evdt_(PRMS.event_prms),
    seed_tracker_(PRMS.seed_prms),
    state_(State::INACTIVE) {

    load_static();

    PathBuffer::MAX_PATH_LEN = PRMS.seed_len;

    for (u64 t = 0; t < EventType::NUM_TYPES; t++) {
        PathBuffer::TYPE_ADDS[t] = t << ((PathBuffer::MAX_PATH_LEN-2)*TYPE_BITS);
    }
    PathBuffer::TYPE_MASK = (u8) ((1 << TYPE_BITS) - 1);

    kmer_probs_ = std::vector<float>(kmer_count<KLEN>());
    prev_paths_ = std::vector<PathBuffer>(PRMS.max_paths);
    next_paths_ = std::vector<PathBuffer>(PRMS.max_paths);
    sources_added_ = std::vector<bool>(kmer_count<KLEN>(), false);

    prev_size_ = 0;
    event_i_ = 0;
    seed_tracker_.reset();

    norm_.set_target(model.get_means_mean(), model.get_means_stdv());
}

Mapper::Mapper(const Mapper &m) : Mapper() {}

Mapper::~Mapper() {

    #ifdef DEBUG_SEEDS
    seeds_out_.close();
    #endif

    for (u32 i = 0; i < next_paths_.size(); i++) {
        next_paths_[i].free_buffers();
        prev_paths_[i].free_buffers();
    }
}

//TODO load from conf?
//auto load in constructor?
void Mapper::load_static() {

    if (fmi.is_loaded()) return;

    fmi.load_index(PRMS.bwa_prefix);
    if (!fmi.is_loaded()) {
        std::cerr << "Error: failed to load BWA index\n";
        abort();
    }

    std::ifstream param_file(PRMS.bwa_prefix + INDEX_SUFF);
    if (!param_file.is_open()) {
        std::cerr << "Error: failed to load uncalled index\n";
        abort();
    }

    std::string param_line;

    char *idx_preset_c = (char *) PRMS.idx_preset.c_str();
    prob_threshes_.resize(64);

    //TODO: clean up parser
    //maybe use toml?
    //try making backwards compatible?
    //definitely more error checking
    while (getline(param_file, param_line)) {
        char *param_name = strtok((char *) param_line.c_str(), "\t");
        char *fn_str = strtok(NULL, "\t");
        //char *path_str = strtok(NULL, "\t");
        if ( !PRMS.idx_preset.empty() && strcmp(param_name, idx_preset_c) ) {
                continue;
        }

        u8 fmbin = prob_threshes_.size() - 1;
        char *prob_str;
        while ( (prob_str = strtok(fn_str, ",")) != NULL ) {
            fn_str = NULL;
            prob_threshes_[fmbin] = atof(prob_str);
            fmbin--;
        }

        for (;fmbin < prob_threshes_.size(); fmbin--) {
            prob_threshes_[fmbin] = prob_threshes_[fmbin+1];
        }
    }

}

float Mapper::get_prob_thresh(u64 fmlen) const {
    return prob_threshes_[__builtin_clzll(fmlen)];
}

float Mapper::get_source_prob() const {
    return prob_threshes_.front();
}

u16 Mapper::get_max_events() const {
    if (event_i_ + PRMS.evt_batch_size > PRMS.max_events) 
        return PRMS.max_events - event_i_;
    return PRMS.evt_batch_size;
}

ReadBuffer &Mapper::get_read() {
    return read_;
}

void Mapper::deactivate() {
    state_ = State::INACTIVE;
    reset_ = false;
}

Paf Mapper::map_read() {
    if (read_.loc_.is_mapped()) return read_.loc_;

    map_timer_.reset();

    norm_.set_signal(evdt_.get_means(read_.full_signal_));

    while (!map_next()) {}

    read_.loc_.set_float(Paf::Tag::MAP_TIME, map_timer_.get());

    return read_.loc_;
}

void Mapper::new_read(ReadBuffer &r) {
    read_.clear();//TODO: probably shouldn't auto erase previous read
    read_.swap(r);
    reset();

    #ifdef DEBUG_SEEDS
    seeds_out_.open(read_.id_ + "_seeds.bed");
    #endif

}

void Mapper::new_read(Chunk &chunk) {
    //if (!chunk_mtx_.try_lock()) {
    //    std::cerr << "Error: failed to create new read '" << read_.id_ << "'\n";
    //    return;
    //}

    if (prev_unfinished(chunk.get_number())) {
        std::cerr << "Error: possibly lost read '" << read_.id_ << "'\n";
    }
    read_ = ReadBuffer(chunk);
    reset();

    //chunk_mtx_.unlock();
}

void Mapper::reset() {
    #ifdef DEBUG_SEEDS
    seeds_out_.close();
    #endif

    prev_size_ = 0;
    event_i_ = 0;
    reset_ = false;
    last_chunk_ = false;
    state_ = State::MAPPING;
    norm_.skip_unread();

    seed_tracker_.reset();
    evdt_.reset();

    chunk_timer_.reset();
    map_timer_.reset();
    map_time_ = 0;
    wait_time_ = 0;
}

u32 Mapper::prev_unfinished(u32 next_number) const {
    return state_ == State::MAPPING && read_.number_ != next_number;
}

bool Mapper::finished() const {
    return state_ == State::SUCCESS || state_ == State::FAILURE;
}

void Mapper::skip_events(u32 n) {
    event_i_ += n;
    prev_size_ = 0;
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
    return read_.chunk_processed_;
}

Mapper::State Mapper::get_state() const {
    return state_;
}

bool Mapper::add_chunk(Chunk &chunk) {
    if (!chunk_mtx_.try_lock()) return false;

    if (!is_chunk_processed() || finished() || reset_) { 
        chunk_mtx_.unlock();
        return false;
    }

    if (read_.chunks_maxed()) {

        set_failed();
        chunk.clear();

        chunk_mtx_.unlock();
        return true;
    }

    bool added = read_.add_chunk(chunk);
    if (added) {
        chunk_timer_.reset();
    }

    chunk_mtx_.unlock();
    return added;
}

u16 Mapper::process_chunk() {
    if (read_.chunk_processed_ || reset_ || 
        !chunk_mtx_.try_lock()) return 0; 

    if (read_.chunk_count() == 1) {
        read_.loc_.set_float(Paf::Tag::QUEUE_TIME, map_timer_.lap());
    }

    wait_time_ += map_timer_.lap();


    float mean;
    //std::cerr << "# got " << read_.get_id() 
    //          << " " << read_.chunk_count() << "\n";

    u16 nevents = 0;
    for (u32 i = 0; i < read_.chunk_.size(); i++) {
        if (evdt_.add_sample(read_.chunk_[i])) {
            mean = evdt_.get_mean();

            if (!norm_.push(mean)) {

                u32 nskip = norm_.skip_unread(nevents);
                skip_events(nskip);

                std::cerr << "#SKIP "
                          << read_.get_id() << " "
                          << nskip << "\n";

                if (!norm_.push(mean)) {
                    map_time_ += map_timer_.lap();

                    chunk_mtx_.unlock();
                    return nevents;
                }
            }
            nevents++;
        }
    }

    read_.chunk_.clear();

    read_.chunk_processed_ = true;

    map_time_ += map_timer_.lap();

    chunk_mtx_.unlock();
    return nevents;
}

void Mapper::set_failed() {
    state_ = State::FAILURE;
    reset_ = false;

    read_.loc_.set_float(Paf::Tag::MAP_TIME, map_time_);
    read_.loc_.set_float(Paf::Tag::WAIT_TIME, wait_time_);
}

bool Mapper::chunk_mapped() {
    //std::cerr << "#check " << read_.get_id() << " " << read_.chunk_processed_ << " " << norm_.empty() << "\n";

    return read_.chunk_processed_ && norm_.empty();
}

bool Mapper::map_chunk() {
    wait_time_ += map_timer_.lap();

    if (reset_ || chunk_timer_.get() > PRMS.chunk_timeout) {
        set_failed();
        read_.loc_.set_ended();
        //std::cerr << "# END timer or reset\n";
        return true;

    } else if (norm_.empty() && 
               read_.chunk_processed_ && 
               read_.chunks_maxed()) {

        chunk_mtx_.lock();

        if (norm_.empty() && read_.chunk_processed_) {
            set_failed();
            chunk_mtx_.unlock();
            return true;
        }

        chunk_mtx_.unlock();
    }

    if (norm_.empty()) {
        //std::cerr << "# stuck empty: " 
        //          << chunk_timer_.get() << " < "
        //          << PRMS.chunk_timeout << "\n";
        return false;
    }


    u16 nevents = get_max_events();
    float tlimit = PRMS.evt_timeout * nevents;

    for (u16 i = 0; i < nevents && !norm_.empty(); i++) {
        if (map_next()) {
            read_.loc_.set_float(Paf::Tag::MAP_TIME, map_time_+map_timer_.get());
            read_.loc_.set_float(Paf::Tag::WAIT_TIME, wait_time_);
            norm_.skip_unread();
            return true;
        }

        if (map_timer_.get() > tlimit) {
            //std::cerr << "#event timeout "
            //          << map_timer_.get() << " "
            //          << tlimit << "\n";
            break;
        }
    }

    map_time_ += map_timer_.lap();

    return false;
}

bool Mapper::map_next() {
    if (norm_.empty() || reset_ || event_i_ >= PRMS.max_events) {
        state_ = State::FAILURE;
        return true;
    }

    float event = norm_.pop();

    //TODO: store kmer_probs_ in static array
    for (u16 kmer = 0; kmer < kmer_probs_.size(); kmer++) {
        kmer_probs_[kmer] = model.match_prob(event, kmer);
    }

    Range prev_range;
    u16 prev_kmer;
    float evpr_thresh;
    bool child_found;

    auto next_path = next_paths_.begin();

    
    //Find neighbors of previous nodes
    for (u32 pi = 0; pi < prev_size_; pi++) {
        if (!prev_paths_[pi].is_valid()) {
            continue;
        }

        child_found = false;

        PathBuffer &prev_path = prev_paths_[pi];
        Range &prev_range = prev_path.fm_range_;
        prev_kmer = prev_path.kmer_;

        evpr_thresh = get_prob_thresh(prev_range.length());
        //evpr_thresh = PRMS.get_path_thresh(prev_path.total_match_len_);

        if (prev_path.consec_stays_ < PRMS.max_consec_stay && 
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

        //Add all the neighbors
        for (u8 b = 0; b < BASE_COUNT; b++) {
            u16 next_kmer = kmer_neighbor<KLEN>(prev_kmer, b);

            if (kmer_probs_[next_kmer] < evpr_thresh) {
                continue;
            }

            Range next_range = fmi.get_neighbor(prev_range, b);

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


        if (!child_found && !prev_path.sa_checked_) {

            //Add seeds for non-extended paths
            //Extended paths will be updated after sources filled in
            update_seeds(prev_path, true);

            #ifdef DEBUG_SEEDS
            print_debug_seeds(prev_path);
            #endif
        }

        if (next_path == next_paths_.end()) {
            break;
        }
    }

    //Create sources between gaps
    if (next_path != next_paths_.begin()) {

        u32 next_size = next_path - next_paths_.begin();

        pdqsort(next_paths_.begin(), next_path);
        //std::sort(next_paths_.begin(), next_path);

        u16 source_kmer;
        prev_kmer = kmer_probs_.size(); 

        Range unchecked_range, source_range;

        for (u32 i = 0; i < next_size; i++) {
            source_kmer = next_paths_[i].kmer_;

            //Add source for beginning of kmer range
            if (source_kmer != prev_kmer &&
                next_path != next_paths_.end() &&
                kmer_probs_[source_kmer] >= get_source_prob()) {

                sources_added_[source_kmer] = true;

                source_range = Range(fmi.get_kmer_range(source_kmer).start_,
                                     next_paths_[i].fm_range_.start_ - 1);

                if (source_range.is_valid()) {
                    next_path->make_source(source_range,
                                           source_kmer,
                                           kmer_probs_[source_kmer]);
                    next_path++;
                }                                    

                unchecked_range = Range(next_paths_[i].fm_range_.end_ + 1,
                                        fmi.get_kmer_range(source_kmer).end_);
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
                kmer_probs_[source_kmer] >= get_source_prob()) {
                
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

            update_seeds(next_paths_[i], false);

            #ifdef DEBUG_SEEDS
            print_debug_seeds(next_paths_[i]);
            #endif
        }
    }

    for (u16 kmer = 0; 
         kmer < kmer_probs_.size() && 
            next_path != next_paths_.end(); 
         kmer++) {

        Range next_range = fmi.get_kmer_range(kmer);

        if (!sources_added_[kmer] && 
            kmer_probs_[kmer] >= get_source_prob() &&
            next_path != next_paths_.end() &&
            next_range.is_valid()) {

            //TODO: don't write to prob buffer here to speed up source loop
            next_path->make_source(next_range, kmer, kmer_probs_[kmer]);
            next_path++;

        } else {
            sources_added_[kmer] = false;
        }
    }

    prev_size_ = next_path - next_paths_.begin();
    prev_paths_.swap(next_paths_);

    //Update event index
    event_i_++;

    SeedGroup sg = seed_tracker_.get_final();

    if (sg.is_valid()) {
        state_ = State::SUCCESS;
        set_ref_loc(sg);

        //#ifdef DEBUG_SEEDS
        //for (u32 pi = 0; pi < prev_size_; pi++) {
        //    print_debug_seeds(prev_paths_[pi]);
        //}
        //#endif

        return true;
    }

    return false;
}

void Mapper::update_seeds(PathBuffer &p, bool path_ended) {

    if (p.is_seed_valid(path_ended)) {

        p.sa_checked_ = true;

        for (u64 s = p.fm_range_.start_; s <= p.fm_range_.end_; s++) {

            //Reverse the reference coords so they both go L->R
            u64 ref_en = fmi.size() - fmi.sa(s) + 1;

            seed_tracker_.add_seed(ref_en, p.match_len(), event_i_ - path_ended);
        }
    }
}

#ifdef DEBUG_SEEDS
void Mapper::print_debug_seeds(PathBuffer &p) {

    if (!p.is_seed_valid(true)) return;

    for (u64 s = p.fm_range_.start_; s <= p.fm_range_.end_; s++) {
        //Reverse the reference coords so they both go L->R
        u64 ref_en = fmi.size() - (fmi.sa(s) + 1);

        bool fwd = ref_en < fmi.size() / 2;

        u64 sa_st;
        if (fwd) sa_st = ref_en - (p.match_len() + KLEN - 1)  + 1;
        else     sa_st = fmi.size() - ref_en - 1;

        std::string rf_name;
        u64 rf_st;
        fmi.translate_loc(sa_st, rf_name, rf_st);

        if (rf_st > fmi.size()) {
            rf_st = 0;
        }
          
        seeds_out_ << rf_name << "\t"
                   << rf_st << "\t"
                   << (rf_st + p.match_len() + KLEN) - 1) << "\t"
                   << event_i_ << "\t"
                   << (fwd ? "+" : "-") << "\n";
    }
}
#endif

u32 Mapper::event_to_bp(u32 evt_i, bool last) const {
    return (evt_i * evdt_.mean_event_len() * ReadBuffer::PRMS.bp_per_samp()) + last*(KLEN - 1);
}                  

void Mapper::set_ref_loc(const SeedGroup &seeds) {
    bool fwd = seeds.ref_st_ < fmi.size() / 2;

    u64 sa_st;
    if (fwd) sa_st = seeds.ref_st_;
    else      sa_st = fmi.size() - (seeds.ref_en_.end_ + KLEN - 1);
    
    std::string rf_name;
    u64 rd_st = event_to_bp(seeds.evt_st_),
        rd_en = event_to_bp(seeds.evt_en_ + PRMS.seed_len, true),
        rf_st = 0,
        rf_len = fmi.translate_loc(sa_st, rf_name, rf_st), //sets rf_st
        rf_en = rf_st + (seeds.ref_en_.end_ - seeds.ref_st_ + KLEN);

    u16 match_count = seeds.total_len_ + KLEN - 1;

    read_.loc_.set_mapped(rd_st, rd_en, rf_name, rf_st, rf_en, rf_len, fwd, match_count);
}


