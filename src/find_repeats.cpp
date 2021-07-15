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

#include <iostream>
#include <string>
#include <list>
#include <cstdlib>
#include "range.hpp"
#include "bwa_fmi.hpp"
#include "self_align_ref.hpp"

int main(int argc, char** argv) {
    std::string bwa_prefix(argv[1]),
                fasta_fname(argv[2]);
    u32 min_k = atoi(argv[3]);

    RefIndex fmi(bwa_prefix);
    

    //for (u64 i = 0; i < 50; i++) {
    //    std::cout << BASE_CHARS[fmi.get_base(i)];
    //}
    //std::cout << "\n";

    //return 0;

    //std::vector< std::vector<u8> > seqs;
    //std::ifstream fasta_in(fasta_fname);
    //std::string fasta_line;
    //while (getline(fasta_in, fasta_line)) {
    //    if (fasta_line[0] == '>') {
    //        seqs.push_back(std::vector<u8>());
    //    } else {
    //        for (char c : fasta_line) {
    //            seqs.back().push_back(BASE_COMP_B[BASE_BYTES[(u8)c]]);
    //        }
    //    }
    //}
    //

    fmi.load_pacseq();
    auto seqs = fmi.get_seqs();
    u64 st = 0;

    for (auto s : seqs) {
        std::string name = s.first;
        u64 len = s.second;

        for (u64 i = 0; i < len; i++) {
            Range r = fmi.get_base_range(fmi.get_base(st+i));
            u64 j = i+1;
            for (; j < len && r.length() > 1; j++) {
                r = fmi.get_neighbor(r, fmi.get_base(st+j));
            }
            j--;

            if (j-i >= min_k) {
                std::cout << (j-i) << "\t" << name << "\t" << i << "\t" << j << "\t";
                for (u32 l = i; l < j; l++) std::cout << BASE_CHARS[BASE_COMP_B[fmi.get_base(st+l)]];
                std::cout << "\n";
            }
        }

        st += len;
    }
}

