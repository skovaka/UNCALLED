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

#ifndef _INCL_KMER_MODEL
#define _INCL_KMER_MODEL

#include <array>
#include <utility>
#include <cmath>
#include "event_detector.hpp"
#include "util.hpp"
#include "bp.hpp"

template<KmerLen KLEN>
class PoreModel {
    public:

    PoreModel() 
        :  loaded_(false) {}

    //TODO: clean up IO
    //maybe load from toml and/or header file
    //make scripts for model to toml and header?
    PoreModel(std::string model_fname, bool complement) {
        kmer_count_ = kmer_count<KLEN>();

        std::ifstream model_in(model_fname);

        //Read header and count number of columns
        std::string header;
        std::getline(model_in, header);
        u8 num_cols = 0;
        bool prev_ws = true;
        for (u8 i = 0; i < header.size(); i++) {
            if (header[i] == ' ' || header[i] == '\t') {
                prev_ws = true;
            } else {
                if (prev_ws)
                    num_cols++;
                prev_ws = false;
            }
        }

        //Variables for reading model
        std::string kmer_str, neighbor_kmer;
        u16 kmer;
        float lv_mean, lv_stdv, sd_mean, sd_stdv, lambda, weight;

        //Check if table includes "ig_lambda" column
        bool has_lambda = num_cols >= 7;

        //Read first line
        if (has_lambda) {
            model_in >> kmer_str >> lv_mean >> lv_stdv >> sd_mean 
                >> sd_stdv >> lambda >> weight;
            lambda_ = lambda;
        } else if (num_cols == 3) {
            model_in >> kmer_str >> lv_mean >> lv_stdv; 
            lambda_ = -1;
        } else if (num_cols == 4) {
            model_in >> kmer_str >> lv_mean >> lv_stdv >> sd_mean; 
            lambda_ = -1;
        } else {
            model_in >> kmer_str >> lv_mean >> lv_stdv >> sd_mean 
                >> sd_stdv >> weight;
            lambda_ = -1;
        }

        //Compute number of kmers (4^k) and reserve space for model


        lv_means_.resize(kmer_count_+1);
        lv_vars_x2_.resize(kmer_count_+1);
        lognorm_denoms_.resize(kmer_count_+1);

        lv_means_[kmer_count_] = 
            lv_vars_x2_[kmer_count_] = 
            lognorm_denoms_[kmer_count_] = -1;

        model_mean_ = 0;

        //Read and store rest of the model
        do {
            //Get unique ID for the kmer
            kmer = str_to_kmer<KLEN>(kmer_str);

            //Complement kmer if needed
            if (complement) {
                kmer = kmer_comp<KLEN>(kmer);
            }

            //Check if kmer is valid
            if (kmer < 0 || kmer >= kmer_count_) {
                std::cerr << "Error: kmer '" << kmer << "' is invalid\n";

                //Store kmer information
            } else {

                //Store model information
                lv_means_[kmer] = lv_mean;
                lv_vars_x2_[kmer] = 2*lv_stdv*lv_stdv;
                lognorm_denoms_[kmer] = log(sqrt(M_PI * lv_vars_x2_[kmer]));

                model_mean_ += lv_mean;
            }

            //Read first line
            if (has_lambda) {
                model_in >> kmer_str >> lv_mean >> lv_stdv >> sd_mean 
                    >> sd_stdv >> lambda >> weight;
                lambda_ = lambda;
            } else if (num_cols == 3) {
                model_in >> kmer_str >> lv_mean >> lv_stdv; 
                lambda_ = -1;
            } else if (num_cols == 4) {
                model_in >> kmer_str >> lv_mean >> lv_stdv >> sd_mean; 
                lambda_ = -1;
            } else {
                model_in >> kmer_str >> lv_mean >> lv_stdv >> sd_mean 
                    >> sd_stdv >> weight;
                lambda_ = -1;
            }

            //Read until eof or correct number of kmers read
        } while (!model_in.eof());

        //Compute model level mean and stdv
        model_mean_ /= kmer_count_;
        model_stdv_ = 0;
        for (u16 kmer = 0; kmer < kmer_count_; kmer++)
            model_stdv_ += pow(lv_means_[kmer] - model_mean_, 2);
        model_stdv_ = sqrt(model_stdv_ / kmer_count_);

        loaded_ = true;
    }

    float match_prob(float samp, u16 kmer) const {
        return (-pow(samp - lv_means_[kmer], 2) / lv_vars_x2_[kmer]) - lognorm_denoms_[kmer];
    }

    float match_prob(const Event &evt, u16 kmer) const {
        return match_prob(evt.mean, kmer);
    }


    float get_mean() const {
        return model_mean_;
    }

    float get_stdv() const {
        return model_stdv_;
    }

    bool is_loaded() const {
        return loaded_;
    }

    private:
    std::vector<float> lv_means_, lv_vars_x2_, lognorm_denoms_;
    float lambda_, model_mean_, model_stdv_;
    u16 kmer_count_;
    bool loaded_, complement_;
};

#endif
