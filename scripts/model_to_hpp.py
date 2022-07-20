#!/usr/bin/env python

import sys
import argparse
import pandas as pd
from uncalled import PoreModel

parser = argparse.ArgumentParser(description="Converts a pore model TSV to a C++ header file")
parser.add_argument("tsv", help="Model TSV file")
parser.add_argument("name", help="Model name")
parser.add_argument("--rna", action="store_true", help="Will set reverse=True")
args = parser.parse_args()

model = PoreModel(args.tsv)

header_name = "_INCL_MODEL_" + args.name.upper()

reverse = "true" if args.rna else "false"

print(f"""\
#ifndef {header_name}
#define {header_name}

#include <vector>

const std::vector<float> model_{args.name}_vals = {{""")

for k in model.KMERS:
    print(f"\t{model.means[k]:.6f}, {model.stdvs[k]:.6f}, //{model.kmer_to_str(k)}")

print(f"""\
}};

const PoreModelPreset model_{args.name} {{
    {{"{args.name}", {reverse}, false, {model.K}, {model.get_kmer_shift(model.K)}}}, model_{args.name}_vals
}};

#endif""")

sys.stderr.write(f"""//Add to pore_model.cpp
#include "models/{args.name}.inl"
model_{args.name}.map()\n""")
