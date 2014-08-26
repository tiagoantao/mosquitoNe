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
        high[dist - 1].append(stat[0])
        point[dist - 1].append(flt(stat[1]))
        low[dist - 1].append(stat[2])
        l = f.readline()
    ax.plot([stats.hmean([[y if y > 0 else 100000 for y in x] for x in
                          point])])

plt.ioff()
numCols = len(ncs)
numRows = len(sampleStrats)
for col in range(len(ncs)):
    for row in range(len(sampleStrats)):
        ax = plt.subplot(numRows, numCols, col * numRows + row + 1)
        doPlot(ax, ncs[col], sampleStrats[row], col == 0, row == numRows)
plt.savefig('fig2.png')
