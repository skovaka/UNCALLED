#include <cmath>
#include "normalizer.hpp"

//const NormalizerParams NORMALIZER_PRMS_DEF = {
//    mode : "model_mom",
//    len : 6000,
//    tgt_mean : 90.20827,
//    tgt_stdv : 12.83266
//};

Normalizer::Normalizer(NormalizerParams p) 
    : PRMS(p),
      signal_(p.len), //TODO need to set
      mean_(0),
      varsum_(0),
      n_(0),
      rd_(0),
      wr_(0),
      is_full_(false),
      is_empty_(true) {
}

Normalizer::Normalizer(float tgt_mean, float tgt_stdv) : Normalizer(NORMALIZER_PRMS_DEF) {
    set_target(tgt_mean, tgt_stdv);
}

void Normalizer::set_target(float mean, float stdv) {
    PRMS.tgt_mean = mean;
    PRMS.tgt_stdv = stdv;
}

void Normalizer::set_signal(const std::vector<float> &signal) {
    signal_ = signal;
    n_ = signal_.size();
    rd_ = wr_ = 0;
    is_full_ = true;
    is_empty_ = false;

    mean_ = 0;
    for (float e : signal_) mean_ += e;
    mean_ /= n_;

    varsum_ = 0;
    for (auto e : signal_) varsum_ += pow(e - mean_, 2);
}

bool Normalizer::push(float newevt) {
    if (is_full_) {
        return false;
    }

    double oldevt = signal_[wr_];
    signal_[wr_] = newevt;

    //Based on https://stackoverflow.com/questions/5147378/rolling-variance-algorithm
    if (n_ == signal_.size()) {
        double oldmean = mean_;
        mean_ += (newevt - oldevt) / signal_.size();
        varsum_ += (newevt + oldevt - oldmean - mean_) * (newevt - oldevt);

    //Based on https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford's_Online_algorithm
    } else {
        n_++;
        double dt1 = newevt - mean_;
        mean_ += dt1 / n_;
        double dt2 = newevt - mean_;
        varsum_ += dt1*dt2;
    }

    wr_ = (wr_ + 1) % signal_.size();

    is_empty_ = false;
    is_full_ = wr_ == rd_;

    return true;
}

void Normalizer::set_length(u32 len) {
    if (len != 0 && len != PRMS.len) {
        PRMS.len = len;
        signal_.resize(len);
    }
}

void Normalizer::reset(u32 buffer_size) {
    n_ = 0;
    rd_ = 0;
    wr_ = 0;
    mean_ = varsum_ = 0;
    is_full_ = false;
    is_empty_ = true;

    set_length(buffer_size);

    signal_[0] = 0;
}

float Normalizer::kmer_current() const {
    return mean_;
}

float Normalizer::kmer_stdv() const {
    return sqrt(varsum_ / n_);
}

float Normalizer::get_scale() const {
    return PRMS.tgt_stdv / kmer_stdv();
}

float Normalizer::get_shift() const {
    return get_shift(get_scale());
}

float Normalizer::get_shift(float scale) const {
    return PRMS.tgt_mean - scale * mean_;
}

float Normalizer::at(u32 i) const {
    float scale = PRMS.tgt_stdv / sqrt(varsum_ / n_);
    float shift = PRMS.tgt_mean - scale * mean_;
    return scale * signal_[i] + shift;
}

float Normalizer::pop() {
    float e = at(rd_);

    rd_ = (rd_+1) % signal_.size();
    is_empty_ = rd_ == wr_;
    is_full_ = false;

    return e;
}

//TODO use mod instead?
u32 Normalizer::unread_size() const {
    if (rd_ < wr_) return wr_ - rd_;
    else return (n_ - rd_) + wr_;
}

u32 Normalizer::skip_unread(u32 nkeep) {
    if (nkeep >= unread_size()) return 0;

    is_full_ = false;
    is_empty_ = nkeep == 0;

    u32 new_rd;
    if (nkeep <= wr_) new_rd = wr_ - nkeep;
    else new_rd = n_ - (nkeep - wr_);

    u32 nskip;
    if (new_rd > rd_) nskip = new_rd - rd_;
    else nskip = (n_ - rd_) + new_rd;

    rd_ = new_rd;
    return nskip;
}

bool Normalizer::empty() const {
    return is_empty_;
}

bool Normalizer::full() const {
    return is_full_;
}
