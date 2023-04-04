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

#ifndef _INCL_CONF
#define _INCL_CONF

#include <iostream>
#include <vector>
#include <cfloat>
#include "dtw.hpp"
#include "toplevel_prms.hpp"

#ifdef PYBIND
#include "pybind11/pybind11.h"
#endif

class Config {

    public:

    bool cprof;

    //ReadBuffer::Params    read_buffer    = ReadBuffer::PRMS;
    NormalizerParams    normalizer     = NORMALIZER_PRMS_DEF;
    EventDetector::Params event_detector = EventDetector::PRMS_DEF;
    PoreModelParams       pore_model     = PORE_MODEL_PRMS_DEF;

    DtwParams             dtw            = DTW_PRMS_DEF;

    Config() : cprof(false) {}

    void set_r94_rna() {
        pore_model.sample_rate = 3012;
        pore_model.bases_per_sec  = 70;
        pore_model.reverse    = true;

        event_detector.window_length1 = 7;
        event_detector.window_length2 = 12;
        event_detector.threshold1     = 2.8;
        event_detector.threshold2     = 18.0;
    }

    //void export_static() {
    //    ReadBuffer::PRMS = read_buffer;
    //}

    EventDetector::Params &get_event_detector_prms() {
        return event_detector;
    }

    void set_event_detector(const EventDetector::Params &p) {
        event_detector = p;
    }

    #ifdef PYBIND
    #define DEF_METH(P, D) c.def(#P, &Config::P, D);
    #define DEFPRP(P) c.def_property(#P, &Config::get_##P, &Config::set_##P);
    #define DEFPRP_DOC(P) c.def_property( \
                #P, \
                &Config::get_##P, \
                &Config::set_##P, \
                Config::doc_##P());

        
    static std::vector<std::string> _PARAM_GROUPS, _GLOBAL_PARAMS;

    static void pybind_defs(pybind11::class_<Config> &c) {
        c.def(pybind11::init());

        DEF_METH(set_r94_rna, "Sets parameters for RNA sequencing")
        //DEF_METH(export_static, "Sets static parameters for ReadBuffer and Mapper")

        #define CONF_ATTR(P, D, V) \
            c.def_readwrite(#P, &Config::P, D); \
            V.push_back(std::string(#P));
        #define CONF_GROUP(P, D) CONF_ATTR(P, D, _PARAM_GROUPS)
        #define CONF_GLOBAL(P, D) CONF_ATTR(P, D, _GLOBAL_PARAMS)

        //CONF_GROUP(read_buffer, "ReadBuffer parameters")

        CONF_GROUP(normalizer, "") 
        CONF_GROUP(event_detector, "")
        CONF_GROUP(pore_model, "") 
        CONF_GROUP(dtw, "")

        c.def_readonly_static("_PARAM_GROUPS", &Config::_PARAM_GROUPS);
        c.def_readonly_static("_GLOBAL_PARAMS", &Config::_GLOBAL_PARAMS);
    }
    #endif
};

#endif
