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

LD = {}
coanc = {}
het = {}
ys = []
for t in range(cfg.seasonGen, cfg.gens):
    y = cfg.A * math.cos(2 * math.pi * (t - cfg.seasonGen) / cfg.T) + cfg.B
    LD[t] = []
    coanc[t] = []
    het[t] = []
    ys.append(y)
    #print t, y

for numIndivs, numLoci in cfg.sampleStrats:
    fname = myUtils.getStatName(cfg, numIndivs, numLoci)
    for rec in myUtils.getStat(open(fname)):
        if rec["type"] == "temp":
            print rec.keys(), rec["g1"], rec["g2"], rec["Pollak"]
            continue
        LD[rec["gen"]].append(rec["LD"][0])
        het[rec["gen"]].append(rec["Het"][-1])
        coanc[rec["gen"]].append(rec["coanc"])

xs = LD.keys()
xs.sort()
lds = []
coancs = []
hets = []
for x in xs:
    lds.append(stats.hmean([v if v >= 0 else 100000 for v in LD[x]]))
    coancs.append(stats.hmean([v if v >= 0 else 100000 for v in coanc[x]]))
    hets.append(stats.hmean([v if v >= 0 else 100000 for v in het[x]]))
pylab.plot(xs, ys, label="Nc")
pylab.plot(xs, lds, label="LD")
pylab.plot(xs, coancs, label="Coanc")
pylab.plot(xs, hets, label="Het")
pylab.legend()
pylab.show()
