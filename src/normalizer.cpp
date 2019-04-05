#include "normalizer.hpp"


Normalizer::Normalizer(const KmerModel &model, 
                       u32 buffer_size) 
    : model_(model),
      events_(buffer_size),
      mean_(0),
      varsum_(0),
      n_(0),
      rd_(0),
      wr_(0),
      is_full_(false),
      is_empty_(true) {

}


bool Normalizer::add_event(float newevt) {
    if (is_full_) {
        //std::cout << "# FULL UP\n";
        return false;
    }

    double oldevt = events_[wr_];
    events_[wr_] = newevt;
    //std::cout << wr_ << " " << rd_ << " " << newevt << " " 
    //          << is_empty_ << " " << is_full_ << " push\n";

    //Based on https://stackoverflow.com/questions/5147378/rolling-variance-algorithm
    if (n_ == events_.size()) {
        double oldmean = mean_;
        mean_ += (newevt - oldevt) / events_.size();
        varsum_ += (newevt + oldevt - oldmean - mean_) * (newevt - oldevt);

    //Based on https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford's_Online_algorithm
    } else {
        n_++;
        double dt1 = newevt - mean_;
        mean_ += dt1 / n_;
        double dt2 = newevt - mean_;
        varsum_ += dt1*dt2;
    }


    wr_ = (wr_ + 1) % events_.size();

    is_empty_ = false;
    is_full_ = wr_ == rd_;

    return true;
}

void Normalizer::reset(u32 buffer_size = 0) {
    n_ = 0;
    rd_ = 0;
    wr_ = 0;
    mean_ = varsum_ = 0;
    is_full_ = false;
    is_empty_ = true;

    if (buffer_size != 0 && buffer_size != events_.size()) {
        events_.resize(buffer_size);
    }

    events_[0] = 0;
}

NormParams Normalizer::get_params() const {
    NormParams p;

    p.scale = model_.model_stdv_ / sqrt(varsum_ / n_);
    p.shift = model_.model_mean_ - p.scale * mean_;

    return p;
}

float Normalizer::pop_event() {
    NormParams np = get_params();
    float e = (float) (np.scale * events_[rd_] + np.shift);
    //std::cout << wr_ << " " << rd_ << " " << e << " " 
    //          << is_empty_ << " " << is_full_ << " pop\n";
    rd_ = (rd_+1) % events_.size();
    is_empty_ = rd_ == wr_;
    is_full_ = false;
    return e;
}

u32 Normalizer::unread_size() const {
    if (rd_ < wr_) return wr_ - rd_;
    else return (n_ - rd_) + wr_;
}

u32 Normalizer::skip_unread(u32 nkeep) {
    if (nkeep >= unread_size()) return 0;

    //std::cout << "# skip " 
    //          << nkeep << " "
    //          << rd_ << " "
    //          << wr_ << " "
    //          << n_ << " ";

    is_full_ = false;
    is_empty_ = nkeep == 0;

    u32 new_rd;
    if (nkeep <= wr_) new_rd = wr_ - nkeep;
    else new_rd = n_ - (nkeep - wr_);
    std::cout << new_rd << " "; 

    u32 nskip;
    if (new_rd > rd_) nskip = new_rd - rd_;
    else nskip = (n_ - rd_) + new_rd;
    std::cout << nskip << "\n"; 

    rd_ = new_rd;
    return nskip;
}

bool Normalizer::empty() const {
    return is_empty_;
}

bool Normalizer::full() const {
    return is_full_;
}
