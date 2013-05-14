import math
import sys

import pylab
from scipy import stats

import myUtils

if len(sys.argv) != 3:
    print "Syntax:", sys.argv[0], "<conffile> <tempStat>"
    sys.exit(-1)

etc = myUtils.getEtc()

cfg = myUtils.getConfig(sys.argv[1])
tempStat = sys.argv[2]

res = {}
ys = []
for t in cfg.futureGens:
    y = cfg.A * math.cos(2 * math.pi * (t - cfg.seasonGen) / cfg.T) + cfg.B
    res[t] = {}
    #print t, y

for numIndivs, numLoci in cfg.sampleStrats:
    fname = myUtils.getStatName(cfg, numIndivs, numLoci)
    for rec in myUtils.getStat(open(fname)):
        if rec["type"] != "temp":
            continue
        g1l = res[rec["g2"]].setdefault(rec["g1"], [])
        g1l.append(rec[tempStat][-1])

xs = res.keys()
xs.sort()
plt = {}
for x in xs:
    for g1 in res[x].keys():
        g1l = plt.setdefault(g1, [])
        g1l.append(stats.hmean([v if v >= 0 else 100000 for v in res[x][g1]]))
print len(xs), plt.keys(), len(plt[56])
pylab.title(tempStat)
for g1, vals in plt.items():
    print g1
    pylab.plot(xs, vals, label=str(g1))
pylab.legend()
pylab.show()
