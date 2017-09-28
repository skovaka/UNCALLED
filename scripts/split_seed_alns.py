import sys

infile = open(sys.argv[1])
out_prefix = sys.argv[2]

aln_id = 0

new_file = True

for line in infile:
    if new_file:
        out = open("%s%03d.txt" % (out_prefix, aln_id), 'w')
        new_file = False

    out.write(line)

    if line[0] == "=" and "done" in line:
        aln_id += 1
        new_file = True
        out.close()

out.close()

