import sys

MER_LEN = 6
BASES = ['A', 'C', 'G', 'T']

def i_to_mer(i):
    ret = [0]*MER_LEN
    for n in range(0, MER_LEN):
        ret[MER_LEN-n-1] = BASES[(i >> n*2) & 3]
    return "".join(ret)

for i in sys.argv[1:]:
    print i_to_mer(int(i))
