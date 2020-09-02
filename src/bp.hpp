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

#ifndef _INCL_BP 
#define _INCL_BP 

#include <string>
#include <vector>
#include <cstdint>
#include <cmath>

#define BASE_COUNT 4

const char BASE_CHARS[] {'A', 'C', 'G', 'T', 'N'};

const u8 BASE_BYTES[] 
     {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //0-15 
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //16-31
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //32-47
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //48-63
      4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, //64-79 (A,C,G)
      4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //80-95 (T)
      4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, //96-111 (a,c,g)
      4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};//112-127 (t)

const char BASE_COMP_C[] 
     {'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N', //0-15  ga
      'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N', //16-31 rb
      'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N', //32-47 ag
      'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N', //48-63 e!
      'N','T','N','G','N','N','N','C','N','N','N','N','N','N','N','N', //64-79 (A,C,G)
      'N','N','N','N','A','N','N','N','N','N','N','N','N','N','N','N', //80-95 (T)
      'N','t','N','g','N','N','N','c','N','N','N','N','N','N','N','N', //96-111 (a,c,g)
      'N','N','N','N','a','N','N','N','N','N','N','N','N','N','N','N'};//112-127 (t)

const u8 BASE_COMP_B[] {3, 2, 1, 0};


enum KmerLen {k2=2, k3=3, k4=4, k5=5};

#define KMASK(k) ( (1 << (2*k)) - 1 )

template <KmerLen k>
inline u16 kmer_count() {
    return (u16) pow((u16)BASE_COUNT, (u16)k);
}

template <KmerLen k>
u16 str_to_kmer(std::string kmer, u64 offset=0) {
    u16 index = BASE_BYTES[(u8) kmer[offset]];
    for (u8 i = 1; i < (u8) k; i++) {
        index = (index << 2) | BASE_BYTES[(u8) kmer[offset+i]];
    }
    return index;
}

template <KmerLen k>
u16 kmer_comp(u16 kmer) {
    return kmer ^ KMASK((u8) k);
}

template <KmerLen k>
u16 kmer_revcomp(u16 kmer) {
    u16 r = ~kmer;
    r = ( (r >> 2 & 0x3333) | (r & 0x3333) << 2 );
    r = ( (r >> 4 & 0x0F0F) | (r & 0x0F0F) << 4 );
    r = ( (r >> 8 & 0x00FF) | (r & 0x00FF) << 8 );
    return r >> (2 * (8 - k));
}

template <KmerLen KLEN>
std::vector<u16> kmers_revcomp(const std::vector<u16> &kmers) {
    std::vector<u16> rev;
    rev.reserve(kmers.size());
    for (auto k = kmers.rbegin(); k != kmers.rend(); k++) {
        rev.push_back(kmer_revcomp<KLEN>(*k));
    }
    return rev;
}

template <KmerLen k>
u8 kmer_head(u16 kmer) {
    return (u8) ((kmer >> (2*( (u8)k ) - 2)) & 0x3);
}

template <KmerLen k>
u16 kmer_neighbor(u16 kmer, u8 i) {
    return ((kmer << 2) & KMASK(k)) | i; 
}

template <KmerLen k>
u8 kmer_base(u16 kmer, u8 i) {
    return (u8) ((kmer >> (2 * ((u16)k-i-1))) & 0x3);
}

template <KmerLen k>
std::string kmer_to_str(u16 kmer) {
    std::string s(k, 'N');
    for (u8 i = 0; i < k; i++) {
        s[i] = BASE_CHARS[kmer_base<k>(kmer, i)];
    }
    return s;
}

template <KmerLen KLEN>
std::vector<u16> seq_to_kmers(u8 *seq, u64 st, u64 en) {
    std::vector<u16> ret;

    u64 pst = st >> 2,
        pen = ((en) >> 2)+1;

    u64 i = 0;
    u16 kmer = 0;
    u8 bst = (st&3), ben;

    for (u64 j = pst; j < pen; j++) {
        ben = j == pen-1 ? (en&3) : 4;
        for (u8 k = bst; k < ben; k++) {
            kmer = kmer_neighbor<KLEN>(kmer, (seq[j] >> ((k^3) << 1) ) & 3);
            if (++i >= KLEN) ret.push_back(kmer);
        }
        bst = 0;
    }

    return ret;
}


#endif
