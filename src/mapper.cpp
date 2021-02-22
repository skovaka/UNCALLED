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
#include <exception>
#include "mapper.hpp"
#include "model_r94.inl"

Mapper::Params Mapper::PRMS {
    seed_len        : 22,
    min_rep_len     : 0,
    max_rep_copy    : 50,
    max_paths       : 10000,
    max_consec_stay : 8,
    max_events      : 30000,
    max_stay_frac   : 0.5,
    min_seed_prob   : -3.75,
    evt_batch_size  : 5,
    evt_timeout     : 10.0,
    chunk_timeout   : 4000.0,
    skip_chunks     : 0,
    bwa_prefix      : "",
    idx_preset      : "default",
    model_path      : "",
    seed_prms       : SeedTracker::PRMS_DEF,
    norm_prms       : Normalizer::PRMS_DEF,
    event_prms      : EventDetector::PRMS_DEF,
    evt_prof_prms : EventProfiler::PRMS_DEF

    #ifdef DEBUG_OUT
    , dbg_prefix : "dbg_"
    #endif
};

BwaIndex<KLEN> Mapper::fmi;
std::vector<float> Mapper::prob_threshes_;

PoreModel<KLEN> Mapper::model = pmodel_r94_complement;

const std::array<u8,Mapper::EVENT_TYPES.size()> Mapper::EVENT_TYPES = {
    Mapper::EVENT_STAY,
    Mapper::EVENT_MOVE
};
u32 Mapper::PATH_MASK = 0;
u32 Mapper::PATH_TAIL_MOVE = 0;

Mapper::Mapper() :
    evdt_(PRMS.event_prms),
    evt_prof_(PRMS.evt_prof_prms),
    norm_(PRMS.norm_prms),
    seed_tracker_(PRMS.seed_prms),
    state_(State::INACTIVE) {

    load_static();

    for (u32 i = 0; i < PRMS.seed_len; i++) {
        PATH_MASK |= 1 << i;
    }
    PATH_TAIL_MOVE = 1 << (PRMS.seed_len-1);

    kmer_probs_ = std::vector<float>(kmer_count<KLEN>());

    PathBuffer::reset_count();//TODO is there a better way?
    prev_paths_ = std::vector<PathBuffer>(PRMS.max_paths);

    PathBuffer::reset_count();
    next_paths_ = std::vector<PathBuffer>(PRMS.max_paths);

    sources_added_ = std::vector<bool>(kmer_count<KLEN>(), false);

    prev_size_ = 0;
    event_i_ = 0;
    seed_tracker_.reset();

    norm_.set_target(model.get_means_mean(), model.get_means_stdv());
}

Mapper::Mapper(const Mapper &m) : Mapper() {}

Mapper::~Mapper() {
    dbg_close_all();

    for (u32 i = 0; i < next_paths_.size(); i++) {
        next_paths_[i].free_buffers();
        prev_paths_[i].free_buffers();
    }
}


void Mapper::load_static() {

    if (fmi.is_loaded()) return;

    if (!PRMS.model_path.empty()) {
        model = PoreModel<KLEN>(PRMS.model_path, true);
    }

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

inline u64 Mapper::get_fm_bin(u64 fmlen) {
    return __builtin_clzll(fmlen);
}

float Mapper::get_prob_thresh(u64 fmlen) const {
    return prob_threshes_[get_fm_bin(fmlen)];
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
    dbg_open_all();
}


void Mapper::new_read(Chunk &chunk) {
    if (prev_unfinished(chunk.get_number())) {
        std::cerr << "Error: possibly lost read '" << read_.id_ << "'\n";
    }

    read_ = ReadBuffer(chunk);
    reset();
}

void Mapper::reset() {
    prev_size_ = 0;
    event_i_ = 0;
    reset_ = false;
    last_chunk_ = false;
    state_ = State::MAPPING;
    norm_.skip_unread();
    //norm_.reset();

    seed_tracker_.reset();
    evdt_.reset();
    evt_prof_.reset();

    chunk_timer_.reset();
    map_timer_.reset();
    map_time_ = 0;
    wait_time_ = 0;

    dbg_close_all();

    #ifdef DEBUG_EVENTS
    dbg_events_.clear();
    #endif

    #ifdef DEBUG_CONFIDENCE
    confident_mapped_ = false;
    #endif
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
        dbg_open_all();
        read_.loc_.set_float(Paf::Tag::QUEUE_TIME, map_timer_.lap());
    }

    wait_time_ += map_timer_.lap();

    u16 nevents = 0;
    //if (read_.chunk_count() > PRMS.skip_chunks) {

    for (u32 i = 0; i < read_.chunk_.size(); i++) {
        if (evdt_.add_sample(read_.chunk_[i])) {

            //Add event to profiler
            //Returns true if next event is not masked
            evt_prof_.add_event(evdt_.get_event());
            
            #ifdef DEBUG_EVENTS
            if (evt_prof_.is_full()) {
                dbg_events_.emplace_back(evt_prof_.anno_event());
            }
            #endif

            if (!evt_prof_.event_ready() || 
                read_.chunk_count() <= PRMS.skip_chunks) continue;

            auto evt_mean = evt_prof_.next_mean();

            if (!norm_.push(evt_mean)) {

                u32 nskip = norm_.skip_unread(nevents);
                skip_events(nskip);

                std::cerr << "#SKIP "
                          << read_.get_id() << " "
                          << nskip << "\n";

                if (!norm_.push(evt_mean)) {
                    map_time_ += map_timer_.lap();

                    chunk_mtx_.unlock();
                    return nevents;
                }
            }

            nevents++;
        }
    }
    //}

    dbg_events_out();

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
    return read_.chunk_processed_ && norm_.empty();
}

bool Mapper::map_chunk() {
    wait_time_ += map_timer_.lap();

    if (reset_ || 
        chunk_timer_.get() > PRMS.chunk_timeout ||
        event_i_ >= PRMS.max_events) {

        set_failed();
        read_.loc_.set_ended();
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

        //evpr_thresh = PRMS.get_path_thresh(prev_path.total_move_len_);

        if (prev_path.consec_stays_ < PRMS.max_consec_stay && 
            kmer_probs_[prev_kmer] >= evpr_thresh) {

            next_path->make_child(prev_path, 
                                  prev_range,
                                  prev_kmer, 
                                  kmer_probs_[prev_kmer], 
                                  EVENT_STAY);
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
                                  EVENT_MOVE);

            child_found = true;

            if (++next_path == next_paths_.end()) {
                break;
            }
        }


        if (!child_found && !prev_path.sa_checked_) {

            //Add seeds for non-extended paths
            //Extended paths will be updated after sources filled in
            update_seeds(prev_path, true);

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

    dbg_paths_out();

    SeedCluster sc = seed_tracker_.get_final();

    if (sc.is_valid()) {

        #ifdef DEBUG_CONFIDENCE
        if (!confident_mapped_) {
            read_.loc_.set_int(Paf::Tag::CONFIDENT_EVENT, evt_prof_.mask_idx_map_[event_i_]);
            confident_mapped_ = true;
            #endif

            
            #ifdef DEBUG_SEEDS
            read_.loc_.set_int(Paf::Tag::SEED_CLUSTER, sc.id_);
            #endif


        #ifdef DEBUG_CONFIDENCE
        }
        #else

        set_ref_loc(sc);
        state_ = State::SUCCESS;
        return true;
        #endif
    }

    //dbg_conf_out();

    //Update event index
    event_i_++;

    return false;
}

void Mapper::update_seeds(PathBuffer &path, bool path_ended) {

    if (!path.is_seed_valid(path_ended)) return;

    //TODO: store actual SA coords?
    //avoid checking multiple times!
    path.sa_checked_ = true;

    for (u64 s = path.fm_range_.start_; s <= path.fm_range_.end_; s++) {

        //TODO: store in buffer, replace sa_checked
        //
        //Reverse the reference coords so they both go L->R
        u64 sa_end = fmi.size() - fmi.sa(s);

        u32 ref_len = path.move_count() + KLEN - 1;
        u64 sa_start = sa_end - ref_len + 1;

        //Add seed and store updated seed cluster
        auto clust = seed_tracker_.add_seed(
            sa_end, 
            path.move_count(), 
            event_i_ - path_ended
        );

        #ifdef DEBUG_SEEDS
        dbg_seeds_out(
            path, 
            clust.id_, 
            event_i_ - path_ended, 
            sa_start, 
            ref_len
        );
        #endif
    }
}


u32 Mapper::event_to_bp(u32 evt_i, bool last) const {
    //TODO store bp_per_samp
    return (evt_i * evdt_.mean_event_len() * ReadBuffer::PRMS.bp_per_samp()) + last*(KLEN - 1) + (PRMS.skip_chunks * ReadBuffer::PRMS.chunk_time * ReadBuffer::PRMS.sample_rate * ReadBuffer::PRMS.bp_per_samp());
}                  

void Mapper::set_ref_loc(const SeedCluster &seeds) {
    bool fwd = seeds.ref_st_ < fmi.size() / 2;

    u64 sa_st;
    if (fwd) sa_st = seeds.ref_st_;
    else      sa_st = fmi.size() - (seeds.ref_en_.end_ + KLEN - 1);
    
    std::string rf_name;
    u64 rd_st = event_to_bp(seeds.evt_st_ - PRMS.seed_len),
        rd_en = event_to_bp(seeds.evt_en_, true),
        rd_len = event_to_bp(event_i_, true),
        rf_st = 0,
        rf_len = fmi.translate_loc(sa_st, rf_name, rf_st), //sets rf_st
        rf_en = rf_st + (seeds.ref_en_.end_ - seeds.ref_st_ + KLEN);

    u16 match_count = seeds.total_len_ + KLEN - 1;

    read_.loc_.set_read_len(rd_len);
    read_.loc_.set_mapped(rd_st, rd_en, rf_name, rf_st, rf_en, rf_len, fwd, match_count);

}

#ifdef DEBUG_OUT
u32 Mapper::PathBuffer::count_ = 0;
#endif

Mapper::PathBuffer::PathBuffer()
    : length_(0),
      prob_sums_(new float[PRMS.seed_len+1]) {

    #ifdef DEBUG_OUT
    id_ = count_++;
    #endif
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
    event_moves_ = EVENT_MOVE;
    seed_prob_ = prob;
    fm_range_ = range;
    kmer_ = kmer;
    sa_checked_ = false;


    //path_type_counts_[EVENT_MOVE] = 1;
    //path_type_counts_[EVENT_STAY] = 0;
    total_move_len_ = 1;

    //TODO: don't write this here to speed up source loop
    prob_sums_[0] = 0;
    prob_sums_[1] = prob;

    #ifdef DEBUG_OUT
    parent_ = PRMS.max_paths;
    #endif
}


void Mapper::PathBuffer::make_child(PathBuffer &p, 
                                    Range &range,
                                    u16 kmer, 
                                    float prob, 
                                    u8 move) {

    u8 stay = 1-move;

    length_ = p.length_ + (p.length_ < PRMS.seed_len);
    fm_range_ = range;
    kmer_ = kmer;
    sa_checked_ = p.sa_checked_;
    event_moves_ = ((p.event_moves_ << 1) | move) & PATH_MASK;
    consec_stays_ = (p.consec_stays_ + stay) * stay;

    total_move_len_ = p.total_move_len_ + move;

    if (p.length_ == PRMS.seed_len) {
        std::memcpy(prob_sums_, &(p.prob_sums_[1]), PRMS.seed_len * sizeof(float));
        prob_sums_[PRMS.seed_len] = prob_sums_[PRMS.seed_len-1] + prob;
        seed_prob_ = (prob_sums_[PRMS.seed_len] - prob_sums_[0]) / PRMS.seed_len;
        event_moves_ |= PATH_TAIL_MOVE;

    } else {
        std::memcpy(prob_sums_, p.prob_sums_, length_ * sizeof(float));
        prob_sums_[length_] = prob_sums_[length_-1] + prob;
        seed_prob_ = prob_sums_[length_] / length_;
    }

    #ifdef DEBUG_OUT
    parent_ = p.id_;
    #endif
}

void Mapper::PathBuffer::invalidate() {
    length_ = 0;
}

bool Mapper::PathBuffer::is_valid() const {
    return length_ > 0;
}

u8 Mapper::PathBuffer::stay_count() const {
    return length_ - move_count();
    //return path_type_counts_[EVENT_MOVE];
}

float Mapper::PathBuffer::prob_head() const {
    return prob_sums_[length_] - prob_sums_[length_-1];

}

u8 Mapper::PathBuffer::move_count() const {
    return __builtin_popcount(event_moves_);
    //return path_type_counts_[EVENT_MOVE];
}

u8 Mapper::PathBuffer::type_head() const {
    //return (event_moves_ >> (PRMS.seed_len-2)) & 1;
    return event_moves_ & 1;
}

u8 Mapper::PathBuffer::type_tail() const {
    //return event_moves_ & 1;
    return (event_moves_ >> (PRMS.seed_len-2)) & 1;
}

bool Mapper::PathBuffer::is_seed_valid(bool path_ended) const {

    //All seeds must be same length
    //and have high probability
    return (length_ == PRMS.seed_len &&
            seed_prob_ >= PRMS.min_seed_prob) && (

               //Must be non repetitive,
               //end in a move
               //and not have too many stays
               (fm_range_.length() == 1 &&
                type_head() == EVENT_MOVE &&
                stay_count() <= PRMS.max_stay_frac * PRMS.seed_len) ||

               //Unless path is terminal,
               //not too repetitive,
               //and not too short
               (path_ended &&
                fm_range_.length() <= PRMS.max_rep_copy &&
                move_count() >= PRMS.min_rep_len)
           );
}


bool operator< (const Mapper::PathBuffer &p1, 
                const Mapper::PathBuffer &p2) {
    return p1.fm_range_ < p2.fm_range_ ||
           (p1.fm_range_ == p2.fm_range_ && 
            p1.seed_prob_ < p2.seed_prob_);
}

void Mapper::dbg_open_all() {
    #ifdef DEBUG_OUT
    if (!dbg_opened_) {

        #ifdef DEBUG_SEEDS
        dbg_open(seeds_out_, "_seeds.bed");
        #endif

        #ifdef DEBUG_PATHS
        dbg_open(paths_out_, "_paths.tsv");
        paths_out_ 
            << "id\t"
            << "parent\t"
            << "fm_start\t"
            << "fm_len\t"
            << "kmer\t"
            << "full_len\t"
            << "match_prob\t"
            << "moves\n";
        #endif

        #ifdef DEBUG_EVENTS
        dbg_open(events_out_, "_events.tsv");
        events_out_ 
            << "start\t"
            << "length\t"
            << "mean\t"
            << "stdv\t"
            << "norm_scale\t"
            << "norm_shift\t"
            << "win_mean\t"
            << "win_stdv\t"
            << "win_mask\n";
        #endif

        //#ifdef DEBUG_CONFIDENCE
        //dbg_open(conf_out_, "_conf.tsv");
        //conf_out_ << "top_conf\t"
        //          << "mean_conf\n";
        //#endif

        dbg_opened_ = true;
    }
    #endif
}

#ifdef DEBUG_OUT
void Mapper::dbg_open(std::ofstream &out, const std::string &suffix) {
    if (out.is_open()) {
        out.close();
    }
    std::string fname = PRMS.dbg_prefix + read_.get_id() + suffix;
    out.open(fname);
    if (!out.is_open()) {
        throw std::invalid_argument("failed to open \"" + fname + "\"\n");
    }
}
#endif

void Mapper::dbg_close_all() {
    #ifdef DEBUG_OUT
    if (dbg_opened_) {
        #ifdef DEBUG_SEEDS
        if (seeds_out_.is_open()) seeds_out_.close();
        #endif

        #ifdef DEBUG_PATHS
        if (paths_out_.is_open()) paths_out_.close();
        #endif

        #ifdef DEBUG_EVENTS
        if (events_out_.is_open()) events_out_.close();
        #endif

        ///#ifdef DEBUG_CONFIDENCE
        ///if (conf_out_.is_open()) conf_out_.close();
        ///#endif

        dbg_opened_ = false;
    }
    #endif
}

//void Mapper::dbg_conf_out() {
//    #ifdef DEBUG_CONFIDENCE
//    if (seed_tracker_.empty() || seed_tracker_.get_top_conf() == 0) return;
//    conf_out_ << evt_prof_.mask_idx_map_[event_i_] << "\t"
//              << seed_tracker_.get_best().id_ << "\t"
//              << seed_tracker_.get_top_conf() << "\t"
//              << seed_tracker_.get_mean_conf() << "\n";
//
//    conf_out_.flush();
//    #endif
//}

void Mapper::dbg_events_out() {
    #ifdef DEBUG_EVENTS
    while(!dbg_events_.empty()) {
        auto e = dbg_events_.front();
        //auto evt = std::get<0>(dbg_events_.front());
        //auto mask = std::get<1>(dbg_events_.front());
        events_out_ 
            << e.evt.start << "\t"
            << e.evt.length << "\t"
            << e.evt.mean << "\t"
            << e.evt.stdv << "\t"
            << norm_.get_scale() << "\t"
            << norm_.get_shift() << "\t"
            << e.win_mean << "\t"
            << e.win_stdv << "\t"
            << e.mask << "\n";
        dbg_events_.pop_front();
    }

    events_out_.flush();
    #endif
}

void Mapper::dbg_seeds_out(
        const PathBuffer &path, 
        u32 clust, 
        u32 evt_end,
        u64 sa_start, 
        u32 ref_len) {
    #ifdef DEBUG_SEEDS

    //TODO de-duplicate code
    //should be storing SA coordinate anyway
    
    //TODO clearly deliniate fm_coord, sa_coord(fw/rv), pacseq_coord, ann_coord

    bool fwd = sa_start < (fmi.size() / 2);

    //TODO change sa_ to clarify unstranded
    u32 sa_half;
    if (fwd) {
        sa_half = sa_start;
    } else {
        sa_half = fmi.size() - (sa_start + ref_len - 1);
    }

    std::string rf_name;
    u64 ref_st = 0;
    fmi.translate_loc(sa_half, rf_name, ref_st);

    seeds_out_ << rf_name << "\t"
               << ref_st << "\t"
               << (ref_st + ref_len) << "\t"

               //name field
               << evt_prof_.mask_idx_map_[evt_end] << ":"
               << path.id_ << ":"
               << clust << "\t"

               << (fwd ? "+" : "-") << "\n";

    seeds_out_.flush();
    #endif

}

void Mapper::dbg_paths_out() {
    #ifdef DEBUG_PATHS
    for (u32 i = 0; i < prev_size_; i++) {
        auto &p = prev_paths_[i];

        u32 evt = evt_prof_.mask_idx_map_[event_i_];

        paths_out_ << evt << ":" 
                   << p.id_ << "\t";

        if (p.parent_ < PRMS.max_paths) {
            paths_out_ << evt_prof_.mask_idx_map_[event_i_-1] << ":" 
                       << p.parent_ << "\t";
        } else {
            paths_out_ << evt << ":" 
                       << p.id_ << "\t";
        }

        paths_out_
            << p.fm_range_.start_ << "\t"
            << p.fm_range_.length() << "\t";

        if (p.is_valid()) {
            paths_out_ << kmer_to_str<KLEN>(p.kmer_) << "\t";
        } else {
            paths_out_ << "NNNNN\t"; //TODO store constant 
        }

        paths_out_ 
            << p.total_move_len_ << "\t"
            << p.prob_head() << "\t";


        if (p.is_valid()) {
            for (u32 i = 0; i < p.length_; i++) {
                paths_out_ << ((p.event_moves_ >> i) & 1);
            }
        } else {
            paths_out_ << 0;
        }

        paths_out_ << "\n";
    }
    #endif
}
