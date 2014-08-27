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


def getPlotConfig(fName):
    config = ConfigParser.ConfigParser()
    config.read(fName)
    return (eval(config.get("params", "nc")),
            eval(config.get("params", "span")),
            eval(config.get("params", "sampleStrat")))


ncs, spans, sampleStrats = getPlotConfig(sys.argv[1])


def doPlot(ax, nc, sampleStrat, startCol, endRow):
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
        rep = int(toks[0])
        ref = int(toks[1])
        gen = int(toks[2])
        dist = gen - ref
        if dist > span:
            l = f.readline()
            continue
        stat = toks[6].split('#')  # Pollak 0.01
        high[dist - 1].append(flt(stat[0]))
        point[dist - 1].append(flt(stat[1]))
        low[dist - 1].append(flt(stat[2]))
        l = f.readline()
    ax.plot([stats.hmean([y if y > 0 else 100000 for y in x]) for x in high])
    ax.plot([stats.hmean([y if y > 0 else 100000 for y in x]) for x in point])
    ax.plot([stats.hmean([y if y > 0 else 100000 for y in x]) for x in low])
    ax.set_ylim(0, 2 * nc)
    ax.get_yaxis().set_ticks([nc // 2, nc, 3 * nc // 2, 2 * nc])
    if not startCol:
        ax.set_yticklabels(['', '', '', ''])
    if endRow:
        ax.set_xticklabels([str(sampleStrat) for sampleStrat in sampleStrats])

plt.ioff()
numCols = len(sampleStrats)
numRows = len(ncs)
fig, axs = plt.subplots(numRows, numCols, sharex=True, figsize=(16, 9))
for col in range(len(sampleStrats)):
    for row in range(len(ncs)):
        ax = axs[row, col]
        doPlot(ax, ncs[row], sampleStrats[col], col == 0, row == numRows)
plt.savefig('fig2.png')
