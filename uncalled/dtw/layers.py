from collections import namedtuple
import pandas as pd
import numpy as np

_Layer = namedtuple("_Layer", ["dtype", "label", "fn", "deps"], defaults=[None,None])

#TODO probably move this to AlnTrack
LAYERS = {
    "seq" : {
        "name" : _Layer(str, "Reference Name",
            lambda track: [track.coords.ref_name]*len(track.layers)),
        "seq.pos" : _Layer("Int64", "Reference Coordinate", 
            lambda track: track.coords.pac_to_pos(track.layer_pacs)),
        "strand" : _Layer(str, "Strand",
            lambda track: track.layer_strands),
        "kmer" : _Layer(str, "Kmer",
            lambda track: track.layer_strands),
        "current" : _Layer(np.float32, "Model Current",
            lambda track: track.layer_strands),
    }, "dtw" : {
        "start" : _Layer("Int32", "Sample Start"),
        "length" : _Layer("Int32", "Sample Length"),
        "start_sec" : _Layer(np.float32, "Time (s)"),
        "length_sec" : _Layer(np.float32, "Dwell (s)"),
        "current" : _Layer(np.float32, "Current (pA)"),
        "current_sd" : _Layer(np.float32, "Stdv (pA)"),
        "kmer" : _Layer("UInt32", "Reference k-mer"),
        "events" : _Layer(np.float32, "Event Count"),
        #"kmer_id" : _Layer(np.uint32, "Binary k-mer ID", 
        #                   lambda track: track.
        "events_log2" : _Layer(np.float32, "Event Count (log2)", 
            lambda track: np.log2(track.layers["dtw","events"]),
            [("dtw","events")]),
        "end" : _Layer("Int32", "Sample End",  
            lambda track: track.layers["dtw","start"] + track.layers["dtw","length"],
            [("dtw", "start"), ("dtw", "length")]),
        "middle" : _Layer(np.float32, "Sample Middle",  
            lambda track: track.layers["dtw","start"] + (track.layers["dtw","length"] / 2),
            [("dtw", "start"), ("dtw", "length")]),
        "length_ms" : _Layer(np.float32, "Dwell (ms)",
            lambda track: 1000 * track.layers["dtw","length"] / track.conf.read_buffer.sample_rate,
            [("dtw", "length")]),
        "model" : _Layer(np.float32, "Model Current",
            lambda track: track.model[track.layers["seq","kmer"]],
            [("dtw", "kmer")]),
        "model_diff" : _Layer(np.float32, "Model pA Diff.",
            lambda track: track.layers["dtw","current"].dropna() - track.model[track.layers["seq","kmer"].dropna()],
            [("dtw", "current"),("seq","kmer")]),
        "abs_diff" : _Layer(np.float32, "Abs. Model Diff.",
            lambda track: (track.layers["dtw","current"] - track.model[track.layers["seq","kmer"]]).abs(),
            [("dtw", "current"),("dtw","kmer")]),
        "base" : _Layer(str, "Reference base",
            lambda track: track.model.kmer_base(track.layers["dtw","kmer"], 2)),
    }, "moves" : {
        "start" : _Layer("Int32", "BC Sample Start"),
        "length" : _Layer("Int32", "BC Sample Length"),
        "start_sec" : _Layer(np.float32, "Time (s)"),
        "length_sec" : _Layer(np.float32, "Dwell (s)"),
        "end" : _Layer("Int32", "Sample End",  
            lambda track: track.layers["moves","start"] + track.layers["moves","length"],
            [("moves", "start"), ("moves", "length")]),
        "middle" : _Layer(np.float32, "Sample Middle",  
            lambda track: track.layers["moves","start"] + (track.layers["moves","length"] / 2),
            [("moves", "start"), ("moves", "length")]),
        #"bp" : _Layer("Int32", "Basecaller Base Index"),
        #"error" : _Layer(str, "Basecalled Alignment Error"),
        "indel" : _Layer("Int32", "Basecalled Alignment Indel"),
    }, "band" : {
        "pac_end" : _Layer("Int32", "Mirror Ref. End"),
        "ref_end" : _Layer("Int32", "Mirror Ref. End",
            lambda track: track.coords.pac_to_pos(track.layers["band","pac_end"]),
            [("band", "pac_end")]),
        "sample_start" : _Layer("Int32", "Raw Sample Start"),
        "sample_end" : _Layer("Int32", "Raw Sample End"),
    }, "cmp" : {
        "aln_a" : _Layer("Int32", "Compare alignment V"),
        "aln_b" : _Layer("Int32", "Compare alignment A"),
        "group_b" : _Layer(str, "Compare type"),
        "jaccard" : _Layer(np.float32, "Jaccard Distance", None, 
            [("dtw", "start"), ("dtw", "end"), ("dtw", "length")]),
        "dist" : _Layer(np.float32, "Mean Ref. Distance", None,
            [("dtw", "start"), ("dtw", "end"), ("dtw", "length")]),
    }, "mvcmp" : {
        "aln_a" : _Layer("Int32", "Compare alignment V"),
        "aln_b" : _Layer("Int32", "Compare alignment A"),
        "group_b" : _Layer(str, "Compare type"),
        "jaccard" : _Layer(np.float32, "Jaccard Distance", None, 
            [("dtw", "start"), ("dtw", "end"), ("dtw", "length"), ("moves", "start"), ("moves", "end"), ("moves", "length")]),
        "dist" : _Layer(np.float32, "Mean Ref. Distance", None,                   
            [("dtw", "start"), ("dtw", "end"), ("dtw", "length"), ("moves", "start"), ("moves", "end"), ("moves", "length")]),
    }
}

LAYER_META = pd.concat([
    pd.concat({
        group : pd.DataFrame(layers, index=_Layer._fields).transpose()
    }, names=("group","layer"))  
    for group, layers in LAYERS.items()
])

LAYER_META["base"] = LAYER_META["fn"].isna()

LAYER_DB_GROUPS = ["dtw", "moves", "cmp", "band"]

def parse_layer(layer):
    
    if isinstance(layer, str):
        spl = layer.split(".")
    elif isinstance(layer, tuple):
        spl = layer
    else:
        raise ValueError("Layer must be string or tuple")

    if len(spl) == 2:
        group,layer = spl
    elif len(spl) == 1:
        if layer in LAYER_META.index.get_level_values(0):
            group = layer
            for layer in LAYER_META.loc[group].query("base").index:
                yield (group, layer)
            return
        group = "dtw"
        layer = spl[0]
    else:
        raise ValueError("Invalid layer specifier \"{layer}\", must contain at most one \".\"")

    if not group in LAYERS:
        opts = ",".join(LAYERS.keys())
        raise ValueError("Invalid layer group \"{group}\". Options: {opts}")

    group_layers = LAYERS[group]

    if not layer in group_layers:
        opts = "\", \"".join(group_layers.keys())
        raise ValueError(f"Invalid layer \"{group}.{layer}\". Options: \"{opts}\"")

    yield (group, layer)

def parse_layers(layers, add_deps=True):
    if layers is None:
        return None

    db_layers = list() 
    fn_layers = list() 

    if isinstance(layers, str):
        layers = layers.split(",")

    ret = list()

    parsed = set()

    for layerstr in layers:
        for layer in parse_layer(layerstr):
            if not layer in parsed:
                parsed.add(layer)

                if add_deps:
                    deps = LAYERS[layer[0]][layer[1]].deps
                    if deps is not None:
                        for dep in deps:
                            if not dep in parsed:
                                parsed.add(dep)
                                yield dep

                yield layer
