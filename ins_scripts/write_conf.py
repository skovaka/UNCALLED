import sys
from collections import defaultdict

fasta_in = sys.argv[1]

ins_seqs = defaultdict(list)

for line in open(fasta_in):
    if line[0] != ">": continue

    seq_name = line[1:].strip()

    tabs = seq_name.split(":")
    st,en = map(int, tabs[-1].split("-"))
    ins_name = tabs[0]

    ins_seqs[ins_name].append( (st, seq_name) )

ups = list()
downs = list()

for ins, seqs in ins_seqs.items():
    if len(seqs) != 2:
        if len(seqs) < 2:
            sys.stderr.write("Error: only one sequence for \"%s\"\n", ins)
        else:
            sys.stderr.write("Error: too many sequences for \"%s\"\n", ins)
        sys.exit(1)

    (_,up_name),(_,down_name) = sorted(seqs)

    ups.append(up_name)
    downs.append(down_name)

print("[realtime]")
print("fwd_tgts = [")
for name in ups:
    print("\t\"%s\"," % name)
print("]")

print("rev_tgts = [")
for name in downs:
    print("\t\"%s\"," % name)
print("]")


