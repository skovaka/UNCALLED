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
    bwa_prefix      : "",
    idx_preset      : "default",
    pore_model      : "r94_dna",//_compl",
    seed_tracker    : SeedTracker::PRMS_DEF,
    normalizer      : Normalizer::PRMS_DEF,
    event_detector  : EventDetector::PRMS_DEF,
    event_profiler  : EventProfiler::PRMS_DEF

    #ifdef DEBUG_OUT
    , meta_prefix : "meta_"
    #endif
};


BwaIndex<KLEN> Mapper::fmi;
std::vector<float> Mapper::prob_threshes_;

PoreModel<KLEN> Mapper::model;// = IS_RNA ? pmodel_r94_rna_templ : pmodel_r94_dna_compl;


const std::array<u8,Mapper::EVENT_TYPES.size()> Mapper::EVENT_TYPES = {
    Mapper::EVENT_STAY,
    Mapper::EVENT_MOVE
};
u32 Mapper::PATH_MASK = 0;
u32 Mapper::PATH_TAIL_MOVE = 0;

Mapper::Mapper() :
    evdt_(PRMS.event_detector),
    evt_prof_(PRMS.event_profiler),
    norm_(PRMS.normalizer),
    seed_tracker_(PRMS.seed_tracker),
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

    norm_.set_target(model.model_mean(), model.model_stdv());
}

Mapper::Mapper(const Mapper &m) : Mapper() {}

Mapper::~Mapper() {
    meta_close_all();

    for (u32 i = 0; i < next_paths_.size(); i++) {
        next_paths_[i].free_buffers();
        prev_paths_[i].free_buffers();
    }
}


void Mapper::load_static() {

    if (fmi.bwt_loaded()) return;

    model = PoreModel<KLEN>(PRMS.pore_model, false, ReadBuffer::PRMS.seq_fwd);

    fmi.load_index(PRMS.bwa_prefix);
    if (!fmi.bwt_loaded()) {
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

inline i64 Mapper::get_fm_bin(i64 fmlen) {
    return __builtin_clzll(fmlen);
}

float Mapper::get_prob_thresh(i64 fmlen) const {
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

Paf Mapper::get_paf() const {
    return out_;
}

void Mapper::deactivate() {
    state_ = State::INACTIVE;
    reset_ = false;
}


Paf Mapper::map_read() {
    if (out_.is_mapped()) return out_;

    map_timer_.reset();

    norm_.set_signal(evdt_.get_means(read_.get_signal()));

    while (!map_next()) {}

    out_.set_float(Paf::Tag::MAP_TIME, map_timer_.get());

    return out_;
}

void Mapper::new_read(ReadBuffer &r) {
    if (prev_unfinished(r.get_number())) {
        std::cerr << "Error: possibly lost read '" << read_.id_ << "'\n";
    }

    read_ = r;
    chunk_processed_ = false;
    out_ = Paf(r.get_id(), r.get_channel(), r.get_start());
    
    reset();
    meta_open_all();
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

    meta_close_all();

    #ifdef PYDEBUG
    meta_ = {};
    #endif

    #ifdef DEBUG_EVENTS
    meta_events_.clear();
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
    return chunk_processed_;
}

Mapper::State Mapper::get_state() const {
    return state_;
}

bool Mapper::add_chunk(ReadBuffer &chunk) {
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

    bool added = read_.add_next_chunk(chunk);
    if (added) {
        chunk_processed_ = false;
        chunk_timer_.reset();
    }

    chunk_mtx_.unlock();
    return added;
}

u16 Mapper::process_chunk() {
    if (chunk_processed_ || reset_ || 
        !chunk_mtx_.try_lock()) return 0; 

    if (read_.get_chunk_count() == 1) {
        meta_open_all();
        out_.set_float(Paf::Tag::QUEUE_TIME, map_timer_.lap());
    }

    wait_time_ += map_timer_.lap();

    u16 nevents = 0;
    for (u32 i = 0; i < read_.size(); i++) {
        if (evdt_.add_sample(read_[i])) {

            #ifdef PYDEBUG
            const auto &evt = evdt_.get_event();
            meta_.events.push_back({
                    evt.start, evt.length, 
                    evt.mean, evt.stdv,
                    0,0,0,0
            });
            //meta_.events.push_back(evdt_.get_event());
            #endif

            //Add event to profiler
            //Returns true if next event is not masked
            evt_prof_.add_event(evdt_.get_event());

            #ifdef PYDEBUG
            u32 evt_i = evt_prof_.get_event_count()-1;
            if (evt_prof_.is_full()) {
                meta_.events[evt_i].prof_mean = evt_prof_.get_win_mean();
                meta_.events[evt_i].prof_stdv = evt_prof_.get_win_stdv();
                meta_.events[evt_i].mask = evt_prof_.event_ready();
            }
            #endif

            if (!evt_prof_.event_ready()) continue;

            #ifdef PYDEBUG
            meta_.event_idxs.push_back(evt_i);
            #endif

            auto evt_mean = evt_prof_.next_mean();

            if (!norm_.push(evt_mean)) {

                u32 nskip = norm_.skip_unread(nevents);
                skip_events(nskip);

                if (!norm_.push(evt_mean)) {
                    map_time_ += map_timer_.lap();

                    chunk_mtx_.unlock();
                    return nevents;
                }
            }

            nevents++;
        }
    }

    meta_events_out();

    read_.clear(); //TODO: this seems weird

    chunk_processed_ = true;

    map_time_ += map_timer_.lap();

    chunk_mtx_.unlock();
    return nevents;
}

void Mapper::set_failed() {
    state_ = State::FAILURE;
    reset_ = false;

    out_.set_float(Paf::Tag::MAP_TIME, map_time_);
    out_.set_float(Paf::Tag::WAIT_TIME, wait_time_);
}

bool Mapper::chunk_mapped() {
    return chunk_processed_ && norm_.empty();
}

bool Mapper::map_chunk() {
    wait_time_ += map_timer_.lap();

    if (reset_ || 
        chunk_timer_.get() > PRMS.chunk_timeout ||
        event_i_ >= PRMS.max_events) {

        set_failed();
        out_.set_ended();
        return true;

    } else if (norm_.empty() && 
               chunk_processed_ && 
               read_.chunks_maxed()) {

        chunk_mtx_.lock();

        if (norm_.empty() && chunk_processed_) {
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
            out_.set_float(Paf::Tag::MAP_TIME, map_time_+map_timer_.get());
            out_.set_float(Paf::Tag::WAIT_TIME, wait_time_);
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

    #ifdef PYDEBUG
    meta_.events[meta_.event_idxs[event_i_]].norm_sig = event;
    #endif

    //TODO: store kmer_probs_ in static array
    for (u16 kmer = 0; kmer < kmer_probs_.size(); kmer++) {
        kmer_probs_[kmer] = model.norm_pdf(event, kmer);
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

        //evpr_thresh = PRMS.get_path_thresh(prev_path.total_moves_);

        if (prev_path.consec_stays_ < PRMS.max_consec_stay && 
            kmer_probs_[prev_kmer] >= evpr_thresh) {

            next_path->make_child(prev_path, 
                                  prev_range,
                                  prev_kmer, 
                                  kmer_probs_[prev_kmer], 
                                  EVENT_STAY);
            #ifdef PYDEBUG
            next_path->evt_ = event_i_;
            #endif

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
            #ifdef PYDEBUG
            next_path->evt_ = event_i_;
            #endif

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
                    #ifdef PYDEBUG
                    next_path->evt_ = event_i_;
                    #endif
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
                    #ifdef PYDEBUG
                    next_path->evt_ = event_i_;
                    #endif
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
            #ifdef PYDEBUG
            next_path->evt_ = event_i_;
            #endif
            next_path++;

        } else {
            sources_added_[kmer] = false;
        }
    }

    //#ifdef PYDEBUG
    //meta_.paths.insert(meta_.paths.end(), next_paths_.begin(), next_path);
    //#endif

    prev_size_ = next_path - next_paths_.begin();
    prev_paths_.swap(next_paths_);

    #ifdef PYDEBUG
    for (u32 i = 0; i < prev_size_; i++) {
        auto &p = prev_paths_[i];
		u32 p_id, p_evt;
        if (p.parent_ < PRMS.max_paths) {
			p_id = p.parent_;
            p_evt = evt_prof_.mask_idx_map_[event_i_-1];
        } else {
			p_id = p.id_;
            p_evt = evt_prof_.mask_idx_map_[event_i_];
        }
        meta_.paths.push_back({
            event       : event_i_,
            id          : p.id_,
            parent      : p.parent_ < PRMS.max_paths ? p.parent_ : p.id_,
            fm_start    : p.fm_range_.start_,
            fm_length   : p.fm_range_.length(),
            kmer        : p.kmer_,
            length      : p.length_,
            total_moves : p.total_moves_,
            norm_pdf  : p.prob_head(),
            seed_prob   : p.seed_prob_,
            moves_pac   : p.event_moves_
        });
    }
    #endif

    meta_paths_out();

    SeedCluster sc = seed_tracker_.get_final();

    if (sc.is_valid()) {

        #ifdef PYDEBUG
        if (meta_.conf_evt == 0) {
            meta_.conf_evt = event_i_;
            meta_.conf_clust = sc.id_;
        }
        #else

        set_ref_loc(sc);
        state_ = State::SUCCESS;
        return true;
        #endif
    }

    //meta_conf_out();

    //Update event index
    event_i_++;

    return false;
}

void Mapper::update_seeds(PathBuffer &path, bool path_ended) {

    if (!path.is_seed_valid(path_ended)) return;

    for (auto fm = path.fm_range_.start_; fm <= path.fm_range_.end_; fm++) {

        //TODO: store in buffer, replace sa_checked
        
        //Reverse the reference coords so they both go L->R
        auto pac_end = fmi.fm_to_refmir(fm); //fmi.size() - fmi.sa(s);

        //Add seed and store updated seed cluster
        auto clust = seed_tracker_.add_seed(
            pac_end, 
            path.move_count(), 
            event_i_ - path_ended
        );

        #ifdef PYDEBUG
        u32 ref_len = path.move_count() + KLEN - 1;
        auto pac_start = pac_end - ref_len;// + 1;

        auto loc = fmi.refmir_to_ref_bound(pac_start, pac_end, read_.PRMS.seq_fwd);

        meta_.seeds.push_back({
            ref_id  : loc.ref_id, 
            start   : static_cast<i64>(loc.start), 
            end     : static_cast<i64>(loc.end), 
            fwd     : loc.fwd,
            event   : event_i_-path_ended, 
            path    : path.id_, 
            cluster : clust.id_
        });
        #endif

        #ifdef DEBUG_SEEDS
        meta_seeds_out(
            path, 
            clust.id_, 
            event_i_ - path_ended, 
            sa_start, 
            ref_len
        );
        #endif
    }

    //TODO: store actual SA coords?
    //avoid checking multiple times!
    path.sa_checked_ = true;
}


u32 Mapper::event_to_bp(u32 evt_i, bool last) const {
    //TODO store bp_per_samp
    return (evt_i * evdt_.mean_event_len() * ReadBuffer::PRMS.bp_per_samp()) + last*(KLEN - 1);
}                  

void Mapper::set_ref_loc(const SeedCluster &seeds) {
    auto loc = fmi.refmir_to_ref_bound(seeds.ref_st_, seeds.ref_en_.end_ + KLEN, read_.PRMS.seq_fwd);
    
    auto rd_st = event_to_bp(seeds.evt_st_ - PRMS.seed_len),
        rd_en = event_to_bp(seeds.evt_en_, true),
        rd_len = event_to_bp(event_i_, true);

    u16 match_count = seeds.total_len_ + KLEN - 1;

    out_.set_read_len(rd_len);
    //TODO clean this up
    out_.set_mapped(rd_st, rd_en, loc.ref_name, loc.start, loc.end, loc.ref_len, loc.fwd, match_count);
}

std::vector<bool> Mapper::unpack_moves(u64 moves, u8 length) {
    std::vector<bool> ret(length);
    for (u32 i = 0; i < length; i++) {
        ret[i] = (moves >> i) & 1;
    }
    return ret;
}

#ifdef PYDEBUG
u32 Mapper::PathBuffer::count_ = 0;
#endif

Mapper::PathBuffer::PathBuffer()
    : length_(0),
      prob_sums_(new float[PRMS.seed_len+1]) {

    #ifdef PYDEBUG
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
    total_moves_ = 1;

    //TODO: don't write this here to speed up source loop
    prob_sums_[0] = 0;
    prob_sums_[1] = prob;

    #ifdef PYDEBUG
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

    total_moves_ = p.total_moves_ + move;

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

    #ifdef PYDEBUG
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

std::vector<bool> Mapper::PathBuffer::get_moves() const { 
    std::vector<bool> ret(length_);
    for (u32 i = 0; i < length_; i++) {
        ret[i] = (event_moves_ >> i) & 1;
    }
    return ret;
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

void Mapper::meta_open_all() {
    #ifdef DEBUG_OUT
    if (!meta_opened_) {

        #ifdef DEBUG_SEEDS
        meta_open(seeds_out_, "_seeds.bed");
        #endif

        #ifdef DEBUG_PATHS
        meta_open(paths_out_, "_paths.tsv");
        paths_out_ 
            << "id\t"
            << "parent\t"
            << "fm_start\t"
            << "fm_len\t"
            << "kmer\t"
            << "full_len\t"
            << "norm_pdf\t"
            << "seed_prob\t"
            << "moves\n";
        #endif

        #ifdef DEBUG_EVENTS
        meta_open(events_out_, "_events.tsv");
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
        //meta_open(conf_out_, "_conf.tsv");
        //conf_out_ << "top_conf\t"
        //          << "mean_conf\n";
        //#endif

        meta_opened_ = true;
    }
    #endif
}

#ifdef DEBUG_OUT
void Mapper::meta_open(std::ofstream &out, const std::string &suffix) {
    if (out.is_open()) {
        out.close();
    }
    std::string fname = PRMS.meta_prefix + read_.get_id() + suffix;
    out.open(fname);
    if (!out.is_open()) {
        throw std::invalid_argument("failed to open \"" + fname + "\"\n");
    }
}
#endif

void Mapper::meta_close_all() {
    #ifdef DEBUG_OUT
    if (meta_opened_) {
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

        meta_opened_ = false;
    }
    #endif
}

//void Mapper::meta_conf_out() {
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

void Mapper::meta_events_out() {
    #ifdef DEBUG_EVENTS
    while(!meta_events_.empty()) {
        auto e = meta_events_.front();
        //auto evt = std::get<0>(meta_events_.front());
        //auto mask = std::get<1>(meta_events_.front());
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
        meta_events_.pop_front();
    }

    events_out_.flush();
    #endif
}

void Mapper::meta_seeds_out(
        const PathBuffer &path, 
        u32 clust, 
        u32 evt_end,
        i64 sa_start, 
        u32 ref_len) {
    #ifdef DEBUG_SEEDS

    //TODO de-duplicate code
    //should be storing SA coordinate anyway

    auto loc = fmi.translate_loc(seeds.ref_st_, seeds.ref_en_.end_ + KLEN, read_.PRMS.seq_fwd);
    
    i64 sa_st;

    bool flip = sa_start >= fmi.size() / 2;
    bool fwd = (!flip && read_.PRMS.seq_fwd) || (flip && !read_.PRMS.seq_fwd);

    u32 sa_half;
    if (flip) sa_half = fmi.size() - (sa_start + ref_len - 1);
    else sa_half = sa_start;

    std::string rf_name;
    i64 ref_st = 0;
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

void Mapper::meta_paths_out() {
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
            << p.total_moves_ << "\t"
            << p.prob_head() << "\t"
            << p.seed_prob_ << "\t";


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
