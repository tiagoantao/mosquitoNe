import math
import sys

import pylab
from scipy import stats

import myUtils

if len(sys.argv) != 3:
    print "Syntax:", sys.argv[0], "<conffile> <refGen>"
    sys.exit(-1)

etc = myUtils.getEtc()

cfg = myUtils.getConfig(sys.argv[1])
refGen = int(sys.argv[2])
tempStats = ["Jorde/Ryman", "Pollak", "Nei/Tajima"]

res = {}
ys = []
for t in cfg.futureGens:
    y = cfg.A * math.cos(2 * math.pi * (t - cfg.seasonGen) / cfg.T) + cfg.B
    res[t] = {}
    #print t, y

for numIndivs, numLoci in cfg.sampleStrats:
    fname = myUtils.getStatName(cfg, numIndivs, numLoci)
    for rec in myUtils.getStat(open(fname)):
        if rec["type"] != "temp" or rec["g1"] != refGen:
            continue

        for tempStat in tempStats:
            g1l = res[rec["g2"]].setdefault(tempStat, [])
            g1l.append(rec[tempStat][-1])

xs = res.keys()
xs.sort()
plt = {}
for x in xs:
    for stat in res[x].keys():
        g1l = plt.setdefault(stat, [])
        g1l.append(stats.hmean([v if v >= 0 else 100000 for v in
                                res[x][stat]]))
print plt.keys()
for stat, vals in plt.items():
    print len(xs), len(vals)
    pylab.plot(xs, vals, label=stat)
pylab.legend()
pylab.show()
