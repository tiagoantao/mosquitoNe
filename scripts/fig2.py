import sys
import ConfigParser

import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

import seaborn as sns


import myUtils
from myUtils import flt

if len(sys.argv) not in [2]:
    print "syntax: %s <plotconfig>" % sys.argv[0]
    sys.exit(-1)


def getPlotConfig(fName):
    config = ConfigParser.ConfigParser()
    config.read(fName)
    return (eval(config.get("params", "nc")),
            eval(config.get("params", "span")),
            eval(config.get("params", "sampleStrat")))


ncs, spans, sampleStrats = getPlotConfig(sys.argv[1])


def doPlot(ax, nc, sampleStrat, startCol):
    numIndivs, numLoci = sampleStrat
    span = 24
    high = [[] for x in range(span)]
    point = [[] for x in range(span)]
    low = [[] for x in range(span)]
    mhigh = [[] for x in range(span)]
    mpoint = [[] for x in range(span)]
    mlow = [[] for x in range(span)]
    do_mlne = False
    if (numIndivs, numLoci) == (60, 20) and nc == 1000:
        do_mlne = True

        class Cfg:
            pass
        cfg = Cfg()  # hard-coded :((((
        cfg.demo = 'constant'
        cfg.numIndivs = numIndivs
        cfg.numLoci = numLoci
        cfg.popSize = nc
        cfg.refGens = [20]
        cfg.futureGens = range(21, 50)
        cfg.reps = 100
        for rep, ref, gen, vals in myUtils.getMLNE(cfg, numIndivs, numLoci):
            dist = gen - ref
            if dist < span:
                mlow[dist - 1].append(vals[0])
                mpoint[dist - 1].append(vals[1])
                mhigh[dist - 1].append(vals[2])

    sampRes = []
    cfg = myUtils.getConfig('simple%d' % nc)
    sampRes.append([])
    #we assume that t0 is file(gen[0])-1
    f = open(myUtils.getStatName(cfg, numIndivs, numLoci))
    l = f.readline()
    while l != "":
        if l.find('temp') > -1:
            break
        l = f.readline()
    f.readline()  # pcrits
    l = f.readline()
    while l != "":
        toks = l.rstrip().split(" ")
        #rep = int(toks[0])
        ref = int(toks[1])
        gen = int(toks[2])
        dist = gen - ref
        if dist > span:
            l = f.readline()
            continue
        stat = toks[-3].split('#')  # Pollak 0.02
        high[dist - 1].append(flt(stat[0]))
        point[dist - 1].append(flt(stat[1]))
        low[dist - 1].append(flt(stat[2]))
        l = f.readline()
    ax.plot([None] + [np.median([y if y > 0 and y < 100000 else 100000 for y in x]) for x in high])
    ax.plot([None] + [np.median([y if y > 0 and y < 100000 else 100000 for y in x]) for x in point],
            'k')
    ax.plot([None] + [np.median([y if y > 0 and y < 100000 else 100000 for y in x]) for x in low])
    if do_mlne:
        ax.plot([None] + [np.median([y if y > 0 and y < 100000 else 100000 for y in x])
                          for x in mhigh], 'k-.')
        ax.plot([None] + [np.median([y if y > 0 and y < 100000 else 100000 for y in x])
                          for x in mpoint], '-.')
        ax.plot([None] + [np.median([y if y > 0 and y < 100000 else 100000 for y in x])
                          for x in mlow], '-.')
    ax.set_ylim(0, 3 * nc)
    ax.get_yaxis().set_ticks([nc // 2, nc, 3 * nc // 2, 2 * nc, 3 * nc])
    ax.axhline(nc)
    ax.set_xlim(0, len(point))
    if not startCol:
        ax.set_yticklabels(['', '', '', ''])

plt.ioff()
numCols = len(sampleStrats)
numRows = len(ncs)
fig, axs = plt.subplots(numRows, numCols, sharex=True, figsize=(16, 9))
for col in range(len(sampleStrats)):
    for row in range(len(ncs)):
        ax = axs[row, col]
        doPlot(ax, ncs[row], sampleStrats[col], col == 0)
        if row == 0:
            ymin, ymax = ax.get_ylim()
            xmin, xmax = ax.get_xlim()
            ni, nl = sampleStrats[col]
            ax.text((xmax - xmin) / 2, ymax, 'I = %d, L = %d' % (ni, nl),
                    va='bottom', ha='center')
fig.tight_layout(h_pad=0, w_pad=0.1)
plt.savefig('fig2.png')
plt.savefig('fig2.eps')
