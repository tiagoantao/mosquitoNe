import math
import sys

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from scipy import stats

import myUtils
from myUtils import flt

if len(sys.argv) != 2:
    print "Syntax:", sys.argv[0], "<conffile>"
    sys.exit(-1)

etc = myUtils.getEtc()

cfg = myUtils.getConfig(sys.argv[1])

LD = {}
coanc = {}
het = {}
temp1 = {}
temp2 = {}
mlne1 = {}
mlne2 = {}
ys = []
for t in cfg.futureGens:
    y = cfg.A * math.cos(2 * math.pi * (t - cfg.seasonGen) / cfg.T) + cfg.B
    LD[t] = []
    coanc[t] = []
    het[t] = []
    temp1[t] = []
    temp2[t] = []
    mlne1[t] = []
    mlne2[t] = []
    ys.append(y)
    #print t, y

numIndivs, numLoci = 60, 20
for irep, ref, gen, vals in myUtils.getMLNE(cfg, numIndivs, numLoci):
    if ref == 50:
        mytemp = mlne1
    elif ref == 56:
        mytemp = mlne2
    else:
        continue
    top, point, bot = vals
    mytemp[gen].append(point)

fname = myUtils.getStatName(cfg, numIndivs, numLoci)
f = open(fname)
f.readline()
l = f.readline()
while l != 'temp\n':
    rep, gen = l.rstrip().split(' ')
    rep = int(rep)
    gen = int(gen)
    l = f.readline()
    toks = l.rstrip().split(' ')
    coanc[gen].append(flt(toks[1]))
    l = f.readline()
    toks = l.rstrip().split(' ')
    stat = toks[2].split('#')  # careful with pcrit
    LD[gen].append(flt(stat[1]))
    l = f.readline()
    toks = l.rstrip().split(' ')
    het[gen].append(flt(toks[1]))  # pcrit...
    l = f.readline()
f.readline()  # pcrit spec
l = f.readline()
while l != '':
    toks = l.rstrip().split(' ')
    ref = int(toks[1])
    if ref == 50:
        mytemp = temp1
    elif ref == 56:
        mytemp = temp2
    else:
        l = f.readline()
        continue
    gen = int(toks[2])
    stat = toks[4].split('#')  # careful with pcrit
    mytemp[gen].append(flt(stat[1]))
    l = f.readline()

xs = LD.keys()
xs.sort()
lds = []
coancs = []
hets = []
temp1s = []
temp2s = []
mlne1s = []
mlne2s = []
fig, ax = plt.subplots()
for x in xs:
    lds.append(stats.hmean([v if v >= 0 else 100000 for v in LD[x]]))
    coancs.append(stats.hmean([v if v >= 0 else 100000 for v in coanc[x]]))
    hets.append(stats.hmean([v if v >= 0 else 100000 for v in het[x]]))
    temp1s.append(stats.hmean([v if v >= 0 else 100000 for v in temp1[x]]))
    temp2s.append(stats.hmean([v if v >= 0 else 100000 for v in temp2[x]]))
    mlne1s.append(stats.hmean([v if v >= 0 else 100000 for v in mlne1[x]]))
    mlne2s.append(stats.hmean([v if v >= 0 else 100000 for v in mlne2[x]]))
ax.plot(xs, ys, label="Nc")
ax.plot(xs, lds, label="LD")
ax.plot(xs, coancs, label="Coanc")
ax.plot(xs, hets, label="Het")
ax.plot(xs, temp1s, label="temp1")
ax.plot(xs, temp2s, label="temp2")
ax.plot(xs, mlne1s, label="mlne1s")
ax.plot(xs, mlne2s, label="mlne2s")
ax.legend()
fig.savefig('wave.png')
