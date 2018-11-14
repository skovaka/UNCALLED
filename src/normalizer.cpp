#include "normalizer.hpp"


Normalizer::Normalizer(const KmerModel &model, 
                       const EventParams &event_params, 
                       u32 buffer_size) 
    : model_(model),
      detector_(event_params),
      sum_(buffer_size),
      sumsq_(buffer_size),
      n_(0),
      rd_(0),
      wr_(0) {

    std::cout << model.model_mean_ << "m " << model.model_stdv_ << "s\n";
}

bool Normalizer::add_sample(float s) {
    Event e = detector_.add_sample(s);
    if (e.length == 0) return false;
    
    wr_ = (wr_ + 1) % sum_.size();

    if (wr_ > 0) {
        sum_[wr_] = sum_[wr_-1] + e.mean;
        sumsq_[wr_] = sumsq_[wr_-1] + e.mean*e.mean;
        std::cout << e.mean << " e ";
    } else {
        sum_[wr_] = sum_.back() + e.mean;
        sumsq_[wr_] = sumsq_.back() + e.mean*e.mean;
    }

    n_ += 1;

    return true;
}

float Normalizer::pop_event() {

    //TODO: return negative if past write

    float mean, var;
    if (n_ < sum_.size()) {
        mean = sum_[wr_] / n_;
        var = (sumsq_[wr_] - (mean*mean*n_)) / n_;
    } else {
        u32 st = (wr_ + 1) % sum_.size();
        mean = (sum_[wr_] - sum_[st]) / (sum_.size() - 1);
        var = (sumsq_[wr_] - sumsq_[st] - mean*mean*(sum_.size()-1)) / (sum_.size()-1);
    }

    float scale = model_.model_stdv_ / sqrt(var),
          shift = model_.model_mean_ - scale*mean;

    u32 e = (rd_+1) % sum_.size();

    float ret = (float) (scale * (sum_[e] - sum_[rd_]) + shift);
    rd_ = e;

    return ret;
}
