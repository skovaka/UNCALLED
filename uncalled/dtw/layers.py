from collections import namedtuple
import pandas as pd
import numpy as np

_Layer = namedtuple("_Layer", ["dtype", "label", "fn", "deps"], defaults=[None,None])

#TODO probably move this to AlnTrack
LAYERS = {
    "ref" : {
        "name" : _Layer(str, "Reference Name",
            lambda track: [track.coords.ref_name]*len(track.layers)),
        "coord" : _Layer("Int64", "Reference Coordinate", 
            lambda track: track.coords.pac_to_ref(track.layer_pacs)),
        "strand" : _Layer(str, "Strand",
            lambda track: track.layer_strands),
    }, "dtw" : {
        "current" : _Layer(float, "Current (pA)"),
        "kmer" : _Layer(str, "Reference k-mer"),
        "start" : _Layer("Int32", "Sample Start"),
        "length" : _Layer("Int32", "Sample Length"),
        "events" : _Layer(float, "Event Count"),
        "events_log2" : _Layer(float, "Event Count (log2)", 
            lambda track: np.log2(track.layers["dtw","events"]),
            [("dtw","events")]),
        "end" : _Layer("Int32", "Sample End",  
            lambda track: track.layers["dtw","start"] + track.layers["dtw","length"],
            [("dtw", "start"), ("dtw", "length")]),
        "middle" : _Layer(float, "Sample Middle",  
            lambda track: track.layers["dtw","start"] + (track.layers["dtw","length"] / 2),
            [("dtw", "start"), ("dtw", "length")]),
        "dwell" : _Layer(float, "Dwell (ms/nt)",
            lambda track: 1000 * track.layers["dtw","length"] / track.conf.read_buffer.sample_rate,
            [("dtw", "length")]),
        "model" : _Layer(float, "Model Current",
            lambda track: track.model[track.layers["dtw","kmer"]],
            [("dtw", "kmer")]),
        "model_diff" : _Layer(float, "Model pA Diff.",
            lambda track: track.layers["dtw","current"] - track.model[track.layers["dtw","kmer"]],
            [("dtw", "current"),("dtw","kmer")]),
        "abs_diff" : _Layer(float, "Abs. Model Diff.",
            lambda track: (track.layers["dtw","current"] - track.model[track.layers["dtw","kmer"]]).abs(),
            [("dtw", "current"),("dtw","kmer")]),
        "base" : _Layer(str, "Reference base",
            lambda track: track.model.kmer_base(track.layers["dtw","kmer"], 2)),
    }, "bcaln" : {
        "start" : _Layer("Int32", "BC Sample Start"),
        "length" : _Layer("Int32", "BC Sample Length"),
        "end" : _Layer("Int32", "Sample End",  
            lambda track: track.layers["bcaln","start"] + track.layers["bcaln","length"],
            [("bcaln", "start"), ("bcaln", "length")]),
        "middle" : _Layer(float, "Sample Middle",  
            lambda track: track.layers["bcaln","start"] + (track.layers["bcaln","length"] / 2),
            [("bcaln", "start"), ("bcaln", "length")]),
        #"bp" : _Layer("Int32", "Basecaller Base Index"),
        #"error" : _Layer(str, "Basecalled Alignment Error"),
        "indel" : _Layer("Int32", "Basecalled Alignment Indel"),
    }, "band" : {
        "pac_end" : _Layer("Int32", "Mirror Ref. End"),
        "ref_end" : _Layer("Int32", "Mirror Ref. End",
            lambda track: track.coords.pac_to_ref(track.layers["band","pac_end"]),
            [("band", "pac_end")]),
        "sample_start" : _Layer("Int32", "Raw Sample Start"),
        "sample_end" : _Layer("Int32", "Raw Sample End"),
    }, "cmp" : {
        "aln_a" : _Layer("Int32", "Compare alignment V"),
        "aln_b" : _Layer("Int32", "Compare alignment A"),
        "group_b" : _Layer(str, "Compare type"),
        "jaccard" : _Layer(float, "Jaccard Distance", None, 
            [("dtw", "start"), ("dtw", "end"), ("dtw", "length")]),
        "mean_ref_dist" : _Layer(float, "Mean Ref. Distance", None,
            [("dtw", "start"), ("dtw", "end"), ("dtw", "length")]),
    }, "bc_cmp" : {
        "aln_a" : _Layer("Int32", "Compare alignment V"),
        "aln_b" : _Layer("Int32", "Compare alignment A"),
        "group_b" : _Layer(str, "Compare type"),
        "jaccard" : _Layer(float, "Jaccard Distance", None, 
            [("dtw", "start"), ("dtw", "end"), ("dtw", "length"), ("bcaln", "start"), ("bcaln", "end"), ("bcaln", "length")]),
        "mean_ref_dist" : _Layer(float, "Mean Ref. Distance", None,                   
            [("dtw", "start"), ("dtw", "end"), ("dtw", "length"), ("bcaln", "start"), ("bcaln", "end"), ("bcaln", "length")]),
    }
}

LAYER_META = pd.concat([
    pd.concat({
        group : pd.DataFrame(layers, index=_Layer._fields).transpose()
    }, names=("group","layer"))  
    for group, layers in LAYERS.items()
])

LAYER_META["base"] = LAYER_META["fn"].isna()

LAYER_DB_GROUPS = ["dtw", "bcaln", "cmp", "band"]

def parse_layer(layer):
    
    if isinstance(layer, str):
        spl = layer.split(".")
    elif isinstance(layer, tuple):
        spl = layer
    else:
        raise ValueError("Layer must be string or tuple")

    if len(spl) == 2:
        group,layer = spl
    #TODO allow for full group specification
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
