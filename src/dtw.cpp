#include "dtw.hpp"

const DtwParams
    DTW_PRMS_DEF = {
        DTWSubSeq::NONE, 1, 1, 1, 0.5, 100, 100, 100, 0, "ref_mom", "guided", "abs_diff", "", false,
    }, DTW_PRMS_EVT_GLOB = {
        DTWSubSeq::NONE, 2, 1, 100, 100, 100, 0, 0, 0, "", "ref_mom", "abs_diff", "", false,
    };

#ifdef PYBIND
py::array_t<Coord> get_guided_bands(PyArray<i32> &bc_refs, PyArray<i32> &bc_samps, PyArray<i32> &event_samps, size_t nbands, i32 shift) {
    //size_t nbands = (bc_refs[bc_refs.size()-1]-bc_refs[0]+1) + event_samps.size();
    //std::cout << "THERE\n";
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

#define PY_DTW_PARAM(P, D) p.def_readwrite(#P, &DtwParams::P, D);
void pybind_dtw(py::module_ &m) {
    PYBIND11_NUMPY_DTYPE(Coord, qry, ref);
    
    DtwDF::pybind<DtwDF>(m, "_DtwDF");

    py::class_<DtwParams> p(m, "DtwParams");
    PY_DTW_PARAM(band_mode, "DTW band mode (\"guided\", \"static\", or \"\"/\"none\")");
    PY_DTW_PARAM(norm_mode, "Normalization method");
    PY_DTW_PARAM(mm2_paf, "Path to minimap2 alignments of basecalled reads in PAF format. Used to determine where each should be aligned. Should include cigar string.");
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
    PyArray<Coord>::pybind(m, "PyArrayCoord");  

    m.attr("DTW_PRMS_DEF") = py::cast(DTW_PRMS_DEF);
    m.attr("DTW_PRMS_EVT_GLOB") = py::cast(DTW_PRMS_EVT_GLOB);
}
#endif
