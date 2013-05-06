import math
import sys

import pylab
from scipy import stats

import myUtils

if len(sys.argv) != 2:
    print "Syntax:", sys.argv[0], "<conffile>"
    sys.exit(-1)

etc = myUtils.getEtc()

cfg = myUtils.getConfig(sys.argv[1])

gens = {}
for t in range(cfg.seasonGen, cfg.gens):
    y = cfg.A * math.cos(2 * math.pi * (t - cfg.seasonGen) / cfg.T) + cfg.B
    gens[t] = []
    #print t, y

for numIndivs, numLoci in cfg.sampleStrats:
    fname = myUtils.getStatName(cfg, numIndivs, numLoci)
    for rec in myUtils.getStat(open(fname)):
        if rec["type"] == "temp":
            continue
        gens[rec["gen"]].append(rec["LD"][0])

xs = gens.keys()
xs.sort()
ys = []
for x in xs:
    ys.append(stats.hmean([v if v >= 0 else 100000 for v in gens[x]]))
pylab.plot(xs, ys)
pylab.show()
