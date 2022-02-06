#include "dtw.hpp"

const DtwParams
    DTW_PRMS_DEF = {
        DTWSubSeq::NONE, 1, 1, 1, 0.5, 100, "guided", "abs_diff", ""
    }, DTW_PRMS_EVT_GLOB = {
        DTWSubSeq::NONE, 2, 1, 100, 0, 0, "", "abs_diff", "",
    };

#ifdef PYBIND
#define PY_DTW_PARAM(P, D) p.def_readwrite(#P, &DtwParams::P, D);
void pybind_dtw(py::module_ &m) {
    py::class_<DtwParams> p(m, "DtwParams");
    PY_DTW_PARAM(band_mode, "DTW band mode (\"guided\", \"static\", or \"\"/\"none\")");
    PY_DTW_PARAM(mm2_paf, "Path to minimap2 alignments of basecalled reads in PAF format. Used to determine where each should be aligned. Should include cigar string.");
    PY_DTW_PARAM(cost_fn, "DTW cost function");
    PY_DTW_PARAM(move_cost, "DTW event move (diagonal) penalty");
    PY_DTW_PARAM(stay_cost, "DTW event stay (horizontal) penalty");
    PY_DTW_PARAM(skip_cost, "DTW event skip (vertical) penalty");
    PY_DTW_PARAM(band_width, "DTW band width");
    PY_DTW_PARAM(band_shift, "DTW band shift");

    m.attr("DTW_PRMS_DEF") = py::cast(DTW_PRMS_DEF);
    m.attr("DTW_PRMS_EVT_GLOB") = py::cast(DTW_PRMS_EVT_GLOB);
}
#endif
