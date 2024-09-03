import liftover
from liftover import get_lifter
from collections import defaultdict
import pybedtools as pbt

def liftover_bedtool(bedtool, gen_from='mm10', gen_to='mm39'):
    converter = get_lifter(gen_from, gen_to)
    ls = []
    c = 0
    for c_og, i in enumerate(bedtool):
        s = converter[i[0]][int(i[1])]
        e = converter[i[0]][int(i[2])]
        if (len(e) == 0) or (len(s) == 0):
            continue
        chrom, s, e = s[0][0], s[0][1], e[0][1]
        s, e = sorted([s, e])
        ls.append([chrom, s, e, c, c_og])
        c += 1
    return pbt.BedTool(ls)

