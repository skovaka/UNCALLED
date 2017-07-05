import sys

model_in = open(sys.argv[1])

events = dict()

for line in model_in:
    if line[0] == "#":
        continue
    tabs = line.strip().split()

    events[tabs[0]] = (tabs[1], tabs[2])

for line in sys.stdin:
    if line[0] == ">":
        continue

    seq = line.strip()

    for i in range(0, len(seq) - 6 + 1):
        print i, seq[i:i+6], "\t".join(events[seq[i:i+6]])
