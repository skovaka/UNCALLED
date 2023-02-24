#include "dtw.hpp"

const DtwParams
    DTW_PRMS_DEF = {
        DTWSubSeq::NONE, 1, 1, 1, 0.5, 100, 100, 100, 1, "ref_mom", "guided", "abs_diff", false,
    }, DTW_PRMS_EVT_GLOB = {
        DTWSubSeq::NONE, 2, 1, 100, 100, 100, 0, 0, 1, "", "ref_mom", "abs_diff", false,
    };

#ifdef PYBIND
//py::array_t<Coord> get_guided_bands(PyArray<i32> &bc_refs, PyArray<i32> &bc_samps, PyArray<i32> &event_samps, size_t nbands, i32 shift) {
py::array_t<Coord> get_guided_bands(py::array_t<i32> &bc_refs_, py::array_t<i32> &bc_samps_, py::array_t<i32> &event_samps_, size_t nbands, i32 shift) {
    //size_t nbands = (bc_refs[bc_refs.size()-1]-bc_refs[0]+1) + event_samps.size();
    PyArray<i32> bc_refs(bc_refs_),
                 bc_samps(bc_samps_), 
                 event_samps(event_samps_);

    auto ret = py::array_t<Coord>(nbands);
    auto bands = PyArray<Coord>(ret);

    size_t i = 0, j = 0, b = 0;
    auto r = bc_refs[0], ref_start = r;
    for (size_t b = 0; b < nbands; b++) {
        bands[b].qry = j + shift;
        bands[b].ref = r - ref_start - shift;

		if (i == bc_refs.size()) j++;
        else if (r < bc_refs[i]) r++;
		else if (j == event_samps.size()-1 || bc_samps[i] < event_samps[j]) {
            i++;
            r++;
        } else j++;
	}

    return ret;
}

constexpr DtwDF::NameArray DtwDF::names;


AlnDF moves_to_aln(py::array_t<bool> moves_py, i32 start, i32 stride) {
    PyArray<bool> moves(moves_py);
    IntervalIndex<i64> index; //({{0, moves.size()}});
    IntervalIndex<i32> samples;

    size_t idx_start = 1, idx_end=1;

    auto end = start + stride;
    for (size_t i = 1; i < moves.size(); i++) {
        if (moves[i]) {
            samples.append(start, end);
            idx_end++;
            start = end;
        }
        end += stride;
    }
    samples.append(start, end);
    index.append(idx_start, idx_end+1);

    return AlnDF(index, samples, {}, {});
}

AlnDF read_to_ref_moves(const AlnDF &read_moves, py::array_t<i64> refs_py, py::array_t<i64> qrys_py, i32 del_max, bool backfill_na=true) {
    PyArray<i64> refs(refs_py), qrys(qrys_py);
    if (refs.size() != qrys.size() || refs.size() == 0) {
        throw std::length_error("reference and query coordinates must be the same non-zero length");
    } 

    IntervalIndex<i64> ref_index; //({{0, moves.size()}});
    IntervalIndex<i32> samples;

    Interval<i64> ref_seg;
    Interval<i32> samps;

    auto ref = refs[0], 
         prev_ref = ref-1;

    ref_seg.start = ref;
    //indel = 0;

    for (size_t i = 0; i < qrys.size(); i++) {
        ref = refs[i];
        
        auto j = read_moves.index.get_index(qrys[i]);
        if (j >= read_moves.samples.interval_count()) {
            continue;
        }
        auto gap = ref - prev_ref;
        if (gap > del_max) {
            ref_seg.end = prev_ref+1;
            ref_index.append(ref_seg);
            ref_seg.start = prev_ref = ref;
            ref_seg.end = ref_seg.NA;
        } else if (gap > 1) {
            auto fill = backfill_na ? samps : Interval<i32>();
            for (auto k = prev_ref+1; k < ref; k++) {
                samples.append(samps);
            }
        }
        prev_ref = ref;
        samps = read_moves.samples.coords[j];
        samples.append(samps);
    }

    for (auto k = prev_ref; k < ref; k++) {
        samples.append(Interval<i32>());
    }

    ref_seg.end = refs[refs.size()-1]+1;
    ref_index.append(ref_seg);
    //samples.append(samps);

    return AlnDF(ref_index, samples, {}, {});
}

#define PY_DTW_PARAM(P, D) p.def_readwrite(#P, &DtwParams::P, D);
void pybind_dtw(py::module_ &m) {
    PYBIND11_NUMPY_DTYPE(Coord, qry, ref);
    
    DtwDF::pybind(m);//<DtwDF>, "_DtwDF");

    py::class_<DtwParams> p(m, "DtwParams");
    PY_DTW_PARAM(band_mode, "DTW band mode (\"guided\", \"static\", or \"\"/\"none\")");
    PY_DTW_PARAM(norm_mode, "Normalization method");
    PY_DTW_PARAM(cost_fn, "DTW cost function");
    PY_DTW_PARAM(move_cost, "DTW event move (diagonal) penalty");
    PY_DTW_PARAM(stay_cost, "DTW event stay (horizontal) penalty");
    PY_DTW_PARAM(skip_cost, "DTW event skip (vertical) penalty");
    PY_DTW_PARAM(del_max, "Will remove reference positions overlapping deletions longer than this");
    PY_DTW_PARAM(ins_max, "Will remove events overlapping insertions longer than this");
    PY_DTW_PARAM(band_width, "DTW band width");
    PY_DTW_PARAM(iterations, "Number of DTW iterations to perform");
    PY_DTW_PARAM(band_shift, "DTW band shift");
    PY_DTW_PARAM(save_bands, "Save DTW band coordinates to database");

    m.def("get_guided_bands", &get_guided_bands);
    m.def("moves_to_aln", &moves_to_aln);
    m.def("read_to_ref_moves", &read_to_ref_moves);
    PyArray<Coord>::pybind(m, "PyArrayCoord");  

    m.attr("DTW_PRMS_DEF") = py::cast(DTW_PRMS_DEF);
    m.attr("DTW_PRMS_EVT_GLOB") = py::cast(DTW_PRMS_EVT_GLOB);
}
#endif
