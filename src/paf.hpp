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

#ifndef _INCL_PAF
#define _INCL_PAF

#include <string>
#include <vector>
#include <iostream>
#include "util.hpp"

class Paf {
    public:

    enum Tag {
        MAP_TIME, 
        WAIT_TIME, 
        QUEUE_TIME, 
        RECEIVE_TIME,
        CHANNEL, 
        EJECT, 
        READ_START, 
        IN_SCAN, 
        TOP_RATIO, 
        MEAN_RATIO,
        ENDED,
        KEEP,
        DELAY,
        SEED_CLUSTER,
        CONFIDENT_EVENT
    };

    Paf();
    Paf(const std::string &rd_name, u16 channel = 0, u64 start_sample = 0);

    bool is_mapped() const;
    bool is_ended() const;
    void print_paf() const;
    void set_read_len(u64 rd_len);
    void set_mapped(u64 rd_st, u64 rd_en, 
                    std::string rf_name,
                    u64 rf_st, u64 rf_en, u64 rf_len,
                    bool fwd, u16 matches);
    void set_ended();
    void set_unmapped();

    void set_int(Tag t, int v);
    void set_float(Tag t, float v);
    void set_str(Tag t, std::string v);

    std::string get_rd_name() {
        return rd_name_;
    }

    #ifdef PYBIND
    #define PY_PAF_METH(P) c.def(#P, &Paf::P);
    #define PY_PAF_TAG(P) t.value(#P, Paf::Tag::P);

    static void pybind_defs(pybind11::class_<Paf> &c) {
        c.def(pybind11::init());
        PY_PAF_METH(print_paf);
        PY_PAF_METH(is_mapped);
        PY_PAF_METH(is_ended);
        PY_PAF_METH(set_int);
        PY_PAF_METH(set_float);
        PY_PAF_METH(set_str);

        pybind11::enum_<Paf::Tag> t(c, "Tag");
        PY_PAF_TAG(MAP_TIME);
        PY_PAF_TAG(EJECT);
        PY_PAF_TAG(IN_SCAN);
        PY_PAF_TAG(ENDED);
        PY_PAF_TAG(KEEP);
        PY_PAF_TAG(DELAY);
        t.export_values();
    }

    #endif

    private:
    static const std::string PAF_TAGS[];

    bool is_mapped_, ended_;
    std::string rd_name_, rf_name_;
    u64 rd_st_, rd_en_, rd_len_,
        rf_st_, rf_en_, rf_len_;
    bool fwd_;
    u16 matches_;

    std::vector< std::pair<Tag, int> > int_tags_;
    std::vector< std::pair<Tag, float> > float_tags_;
    std::vector< std::pair<Tag, std::string> > str_tags_;
};

#endif
