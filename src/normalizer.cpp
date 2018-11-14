#include "normalizer.hpp"


Normalizer::Normalizer(const KmerModel &model, 
                       const EventParams &event_params, 
                       u32 buffer_size) 
    : model_(model),
      detector_(event_params),
      events_(buffer_size),
      mean_(0),
      varsum_(0),
      n_(0),
      rd_(0),
      wr_(0) {

}


bool Normalizer::add_sample(float s) {
    Event e = detector_.add_sample(s);
    if (e.length == 0) return false;

    double oldevt = events_[wr_];
    events_[wr_] = e.mean;
    n_ += 1;

    //Based on https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford's_Online_algorithm
    if (n_ <= events_.size()) {
        double dt1 = e.mean - mean_;
        mean_ += dt1 / n_;
        double dt2 = e.mean - mean_;
        varsum_ += dt1*dt2;

    //Based on https://stackoverflow.com/questions/5147378/rolling-variance-algorithm
    } else {
        double oldmean = mean_;
        mean_ += (e.mean - oldevt) / events_.size();
        //varsum_ += (e.mean - oldmean) * (e.mean - mean_) - (oldevt - oldmean) * (oldevt - mean_);
        varsum_ += (e.mean + oldevt - oldmean - mean_) * (e.mean - oldevt);
    }

    wr_ = (wr_ + 1) % events_.size();

    return true;
}

void Normalizer::reset() {
    n_ = 0;
    rd_ = 0;
    wr_ = 0;
    mean_ = varsum_ = 0;
    events_[0] = 0;
    detector_.reset();
}

NormParams Normalizer::get_params() const {
    NormParams p;
    if (n_ <= events_.size()) {
        p.scale = model_.model_stdv_ / sqrt(varsum_ / n_);
    } else {
        p.scale = model_.model_stdv_ / sqrt(varsum_ / events_.size());
    }

    p.shift = model_.model_mean_ - p.scale * mean_;

    return p;
}

float Normalizer::pop_event() {

    //TODO: return negative if past write
    
    NormParams np = get_params();
    float ret = (float) (np.scale * events_[rd_] + np.shift);

    rd_ = (rd_+1) % events_.size();

    return ret;
}
