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
#include "timer.hpp"
#include "self_align_ref.hpp"

int main(int argc, char** argv) {
    std::string bwa_prefix(argv[1]),
                fasta_fname(argv[2]);
    u32 min_k = atoi(argv[3]);

    BwaFMI fmi(bwa_prefix);

    std::vector< std::vector<u8> > seqs;
    std::ifstream fasta_in(fasta_fname);
    std::string fasta_line;
    while (getline(fasta_in, fasta_line)) {
        if (fasta_line[0] == '>') {
            seqs.push_back(std::vector<u8>());
        } else {
            for (char c : fasta_line) {
                seqs.back().push_back(BASE_COMP_B[BASE_BYTES[(u8)c]]);
            }
        }
    }

    for (auto bases : seqs) {
        for (u64 i = 0; i < bases.size(); i++) {
            Range r = fmi.get_full_range(bases[i]);
            u64 j = i+1;
            for (; j < bases.size() && r.length() > 1; j++) {
                r = fmi.get_neighbor(r, bases[j]);
            }
            j--;

            if (j-i >= min_k) {
                std::cout << (j-i) << "\t" << i << "\t" << j << "\t";
                for (u32 l = i; l < j; l++) std::cout << BASE_CHARS[BASE_COMP_B[bases[l]]];
                std::cout << "\n";
            }
        }
    }
}

