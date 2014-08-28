import sys
import ConfigParser

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

import seaborn as sns

from scipy import stats

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
    span = 24
    high = [[] for x in range(span)]
    point = [[] for x in range(span)]
    low = [[] for x in range(span)]
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
    ax.plot([stats.hmean([y if y > 0 else 100000 for y in x]) for x in high])
    ax.plot([stats.hmean([y if y > 0 else 100000 for y in x]) for x in point])
    ax.plot([stats.hmean([y if y > 0 else 100000 for y in x]) for x in low])
    ax.set_ylim(0, 2 * nc)
    ax.get_yaxis().set_ticks([nc // 2, nc, 3 * nc // 2, 2 * nc])
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
