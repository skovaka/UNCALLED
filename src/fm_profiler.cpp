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

#include <string>
#include <iostream>
#include "bwa_fmi.hpp"
#include "params.hpp"
#include "fm_profiler.hpp"

FMProfiler::FMProfiler() {
    range_counts_.resize(PARAMS.fmi.size());
    kmer_counts_.resize(PARAMS.model.kmer_count());
}

void FMProfiler::add_range(Range r) {
    for (u64 i = r.start_; i <= r.end_; i++) {
        range_counts_[i]++;
    }
}

void FMProfiler::add_kmer(u16 k) {
    kmer_counts_[k]++;
}

void FMProfiler::flush_kmers() {
    for (u64 k = 0; k < kmer_counts_.size(); k++) {
        if (kmer_counts_[k] == 0) continue;

        Range r = PARAMS.kmer_fmranges[k];
        for (u64 i = r.start_; i <= r.end_; i++) {
            range_counts_[i] += kmer_counts_[k];
        }

        kmer_counts_[k] = 0;
    }
}

void FMProfiler::combine(const FMProfiler &p) {
    for (u64 i = 0; i < range_counts_.size(); i++) {
        range_counts_[i] += p.range_counts_[i];
    }

    for (u64 i = 0; i < kmer_counts_.size(); i++) {
        kmer_counts_[i] += p.kmer_counts_[i];
    }
}

void FMProfiler::write(const std::string &fname) {
    flush_kmers();
    std::vector<u32> ref_counts(range_counts_.size());

    for (u64 i = 0; i < ref_counts.size(); i++) {
        ref_counts[PARAMS.fmi.sa(i)] = range_counts_[i];
    }

    std::ofstream out(fname);

    u64 i = 0;
    for (auto seq : PARAMS.fmi.get_seqs()) {
        std::string name = seq.first;
        u64 len = seq.second;
        for (u64 j = 0; j < len; j++) {
            out << name << "\t" 
                << j << "\t" 
                << (j+1) << "\t" 
                << ref_counts[i+j] << "\n";
        }
        i += len;
    }
}
