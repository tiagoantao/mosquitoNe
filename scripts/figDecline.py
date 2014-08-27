import sys
import ConfigParser

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from scipy import stats

import myUtils
from myUtils import flt

if len(sys.argv) not in [2]:
    print "syntax: %s <plotconfig>" % sys.argv[0]
    sys.exit(-1)


ncs = [2000]
numIndivs = 160
numLoci = 50


def doPlot(ax, nc):
    sampRes = []
    cfg = myUtils.getConfig('decl-%d-%d' % (nc, nc // 10))
    sampRes.append([])
    f = open(myUtils.getStatName(cfg, numIndivs, numLoci))
    l = f.readline()
    LD = {}
    coanc = {}
    het = {}
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
    f.readline()  # pcrits
    l = f.readline()
    tdecl = {}
    while l != "":
        toks = l.rstrip().split(" ")
        rep = int(toks[0])
        ref = int(toks[1])
        tdecl.setdefault(ref, {})
        gen = int(toks[2])
        tdecl[ref].setdefault(gen, [])
        stat = toks[6].split('#')  # Pollak 0.01
        tdecl[ref][gen].append(flt(stat[1]))
        l = f.readline()
    ax.plot([stats.hmean([y if y > 0 else 100000 for y in x]) for x in point])

plt.ioff()
numCols = len(ncs)
numRows = len(sampleStrats)
for col in range(len(ncs)):
    for row in range(len(sampleStrats)):
        ax = plt.subplot(numRows, numCols, col * numRows + row + 1)
        doPlot(ax, ncs[col], sampleStrats[row], col == 0, row == numRows)
plt.savefig('fig2.png')
