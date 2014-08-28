from __future__ import division
import sys
import ConfigParser

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

import seaborn as sns

import myUtils

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


def doPlot(ax, nc, span, startCol, endRow):
    sampRes = []
    cfg = myUtils.getConfig('simple%d' % nc)
    for numIndivs, numLoci in sampleStrats:
        print numIndivs, numLoci
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
            if gen - ref != span:
                l = f.readline()
                continue
            stat = toks[-3].split('#')  # Pollak 0.02
            val = float(stat[1])
            sampRes[-1].append([val if val > 0 else 100000])
            l = f.readline()
    sns.boxplot(sampRes, notch=0, sym='', ax=ax)
    ax.set_ylim(0, 2 * nc)
    ax.get_yaxis().set_ticks([nc // 2, nc, 3 * nc // 2, 2 * nc])
    if not startCol:
        ax.set_yticklabels(['', '', '', ''])
    if endRow:
        ax.set_xticklabels([str(sampleStrat) for sampleStrat in sampleStrats])
    ax.axhline(nc)

plt.ioff()
numCols = len(spans)
numRows = len(ncs)
fig, axs = plt.subplots(numRows, numCols, sharex=True, figsize=(16, 9))
for col in range(len(spans)):
    for row in range(len(ncs)):
        ax = axs[row, col]
        doPlot(ax, ncs[row], spans[col], col == 0, row == numRows - 1)
        if row == 0:
            ymin, ymax = ax.get_ylim()
            xmin, xmax = ax.get_xlim()
            ax.text((xmax - xmin) / 2, ymax, 'span = %d' % spans[col],
                    va='bottom', ha='center')
fig.tight_layout(h_pad=0, w_pad=0.1)
plt.savefig('fig1.png')
plt.savefig('fig1.eps')
