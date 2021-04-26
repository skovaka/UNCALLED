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

#ifndef _INCL_TOPLEVEL_PRMS 
#define _INCL_TOPLEVEL_PRMS

typedef struct {
    enum class ActiveChs {FULL, EVEN, ODD, NUM};
    ActiveChs active_chs;

    enum class Mode {DEPLETE, ENRICH, SELECTIVE, NUM};
    Mode realtime_mode;

    std::string host;
    u16 port;
    float duration;

    u32 min_active_reads;
    u32 max_active_reads;
} RealtimeParams;

const RealtimeParams REALTIME_PRMS_DEF = {
    active_chs       : RealtimeParams::ActiveChs::FULL,
    realtime_mode    : RealtimeParams::Mode::ENRICH,
    host             : "127.0.0.1",
    port             : 8000,
    duration         : 72,
    min_active_reads : 0,
    max_active_reads : 512
};

typedef struct {
    std::string ctl_seqsum, unc_seqsum, unc_paf;
    float sim_speed, scan_time, scan_intv_time, ej_time;
} SimulatorParams;

const SimulatorParams SIM_PRMS_DEF = {
    ctl_seqsum     : "",  
    unc_seqsum     : "",
    unc_paf        : "",
    sim_speed      : 1.0,
    scan_time      : 10,
    scan_intv_time : 10.0,
    ej_time        : 5400.0,
};

typedef struct {
    u32 min_active_reads;
} MapOrdParams;

const MapOrdParams MAP_ORD_PRMS_DEF = {
    min_active_reads : 0
};

#endif
