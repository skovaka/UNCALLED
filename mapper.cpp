#include "mapper.hpp"

MapperParams::MapperParams(const BwaFMI &fmi,
                     const KmerModel &model,
                     const std::string &probfn_fname,
                     const EventParams &event_params,
                     u32 seed_len, 
                     u32 min_aln_len,
                     u32 min_rep_len, 
                     u32 max_rep_copy, 
                     u32 max_consec_stay,
                     u32 max_paths, 
                     float max_stay_frac,
                     float min_seed_prob, 
                     float min_mean_conf,
                     float min_top_conf)
        : model_(model),
          fmi_(fmi),
          event_params_(event_params),
          seed_len_(seed_len),
          min_rep_len_(min_rep_len),
          max_rep_copy_(max_rep_copy),
          max_paths_(max_paths),
          max_consec_stay_(max_consec_stay),
          min_aln_len_(min_aln_len),
          max_stay_frac_(max_stay_frac),
          min_seed_prob_(min_seed_prob),
          min_mean_conf_(min_mean_conf),
          min_top_conf_(min_top_conf) {
    
    std::ifstream infile(probfn_fname);
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

float MapperParams::get_prob_thresh(u64 fm_length) {
    auto pr = evpr_threshes_.begin();
    for (auto len = evpr_lengths_.begin(); len != evpr_lengths_.end(); len++) {
        if (fm_length > *len) {
            break;
        }
        pr++;
    }

    return *pr;
}

float MapperParams::get_source_prob() {
    return evpr_threshes_.front();
}

u8 Mapper::PathBuffer::MAX_PATH_LEN = 0, 
   Mapper::PathBuffer::TYPE_MASK = 0;

u64 Mapper::PathBuffer::TYPE_ADDS[EventType::NUM_TYPES];

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
    length_ = p.length_ + 1;

    consec_stays_ = p.consec_stays_;
    fm_range_ = range;
    kmer_ = kmer;
    sa_checked_ = p.sa_checked_;

    std::memcpy(path_type_counts_, p.path_type_counts_, EventType::NUM_TYPES * sizeof(u8));

    if (length_ > MAX_PATH_LEN) {
        std::memcpy(prob_sums_, &(p.prob_sums_[1]), MAX_PATH_LEN * sizeof(float));
        prob_sums_[MAX_PATH_LEN] = prob_sums_[MAX_PATH_LEN-1] + prob;
        seed_prob_ = (prob_sums_[MAX_PATH_LEN] - prob_sums_[0]) / MAX_PATH_LEN;
        path_type_counts_[p.type_tail()]--;
    } else {
        std::memcpy(prob_sums_, p.prob_sums_, (length_ + 1) * sizeof(float));
        prob_sums_[length_] = prob_sums_[length_-1] + prob;
        seed_prob_ = (prob_sums_[length_] - prob_sums_[0]) / length_;
    }

    event_types_ = TYPE_ADDS[type] | (p.event_types_ >> TYPE_BITS);

    path_type_counts_[type]++;

    if (type == EventType::STAY) {
        consec_stays_++;
    } else {
        consec_stays_ = 0;
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

Mapper::Mapper(const MapperParams &ap)
    : params_(ap),
      model_(ap.model_),
      fmi_(ap.fmi_),
      event_detector_(EventDetector(ap.event_params_)),
      seed_tracker_(ap.fmi_.size(),
                    ap.min_mean_conf_,
                    ap.min_top_conf_,
                    ap.min_aln_len_,
                    ap.seed_len_)
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

    reset();

}

Mapper::~Mapper() {
    for (u32 i = 0; i < next_paths_.size(); i++) {
        next_paths_[i].free_buffers();
        prev_paths_[i].free_buffers();
    }
    reset();
}

void Mapper::reset() {
    prev_size_ = 0;
    event_i_ = 0;
    seed_tracker_.reset();
    event_detector_.reset();
    #ifdef DEBUG_TIME
    loop1_time_ = fmrs_time_ = fmsa_time_ = sort_time_ = loop2_time_ = fullsource_time_ = 0;
    #endif
}

ReadAln Mapper::add_samples(const std::vector<float> &samples, const std::string &fast5_name) {
    std::vector<Event> events = event_detector_.get_all_events(samples);
    NormParams norm = model_.get_norm_params(events);
    model_.normalize(events, norm);

    ReadAln aln;

    for (u32 e = 0; e < events.size(); e++) {
        aln = add_event(events[e]);
        if (aln.is_valid()) break;
    }

    if (aln.is_valid()) {

        u64 st, en;
        std::string strand;
        if (aln.ref_st_ < fmi_.size() / 2) {
            st = aln.ref_st_;
            en = aln.ref_en_.end_;
            strand = "rev";
        } else {
            st = fmi_.size() - aln.ref_en_.end_ - 1;
            en = fmi_.size() - aln.ref_st_ - 1;
            strand = "fwd";
        }

        i32 rid = bns_pos2rid(fmi_.bns_, st);
        std::string ref_name(fmi_.bns_->anns[rid].name);
        st -= fmi_.bns_->anns[rid].offset;
        en -= fmi_.bns_->anns[rid].offset;

        std::cout << "aligned\t" 
                 << aln.evt_st_ << "-" << aln.evt_en_  << "\t"
                 << ref_name << "\t" << strand << "\t" 
                 << st << "-" << en << "\t"
                 << aln.total_len_ << "\t";
    } else {
        std::cout << "unaligned\t0-" << events.size() << "\tNULL\tNULL\t0-0\t0\t";
    }
   
    std::cout << fast5_name << std::endl;

    return aln;
}


ReadAln Mapper::add_event(const Event &event
                           #ifdef DEBUG_TIME
                           ,std::ostream &time_out
                           #endif
                           #ifdef DEBUG_SEEDS
                           ,std::ostream &seeds_out
                           #endif
                           ) {


    Range prev_range;
    u16 prev_kmer;
    float evpr_thresh;
    bool child_found;
    std::vector<Seed> seeds;

    #ifdef DEBUG_TIME
    Timer timer;
    #endif

    auto next_path = next_paths_.begin();

    for (u16 kmer = 0; kmer < model_.kmer_count(); kmer++) {
        kmer_probs_[kmer] = model_.event_match_prob(event, kmer);
    }

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

        //Add all the neighbors
        for (u8 i = 0; i < ALPH_SIZE; i++) {
            u16 next_kmer = model_.get_neighbor(prev_kmer, i);

            if (kmer_probs_[next_kmer] < evpr_thresh) {
                continue;
            }

            u8 next_base = model_.get_last_base(next_kmer);

            #ifdef DEBUG_TIME
            loop1_time_ += timer.lap();
            #endif

            Range next_range = fmi_.get_neighbor(prev_range, next_base);

            #ifdef DEBUG_TIME
            fmrs_time_ += timer.lap();
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

        if (!child_found && !prev_path.sa_checked_) {
            #ifdef DEBUG_TIME
            loop1_time_ += timer.lap();
            #endif

            update_seeds(prev_path, seeds, true);

            #ifdef DEBUG_TIME
            fmsa_time_ += timer.lap();
            #endif
        }

        if (next_path == next_paths_.end()) {
            break;
        }
    }

    #ifdef DEBUG_TIME
    loop1_time_ += timer.lap();
    #endif

    if (next_path != next_paths_.begin()) {

        u32 next_size = next_path - next_paths_.begin();

        std::sort(next_paths_.begin(), next_path);

        #ifdef DEBUG_TIME
        sort_time_ += timer.lap();
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
            loop2_time_ += timer.lap();
            #endif

            update_seeds(next_paths_[i], seeds, false);

            #ifdef DEBUG_TIME
            fmsa_time_ += timer.lap();
            #endif
        }
    }

    #ifdef DEBUG_TIME
    loop2_time_ += timer.lap();
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
    fullsource_time_ += timer.lap();
    #endif

    prev_size_ = next_path - next_paths_.begin();
    prev_paths_.swap(next_paths_);

    //Update event index
    event_i_++;

    return seed_tracker_.add_seeds(seeds);
}

void Mapper::update_seeds(PathBuffer &p, 
                           std::vector<Seed> &seeds, 
                           bool path_ended) {

    if (p.is_seed_valid(params_, path_ended)) {
        p.sa_checked_ = true;

        Seed seed(event_i_ - path_ended, params_.seed_len_, p.seed_prob_);
        for (u64 s = p.fm_range_.start_; s <= p.fm_range_.end_; s++) {

            //Reverse the reference coords so they both go L->R
            u64 rev_en = fmi_.size() - fmi_.sa(s) + 1;
            seed.set_ref_range(rev_en, p.match_len());
            seeds.push_back(seed);

            #ifdef DEBUG_SEEDS
            seed.print(seeds_out);
            #endif
        }
    }
}


