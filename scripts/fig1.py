import sys
import ConfigParser

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

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
numCols = len(ncs)
numRows = len(spans)


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
            rep = int(toks[0])
            ref = int(toks[1])
            gen = int(toks[2])
            if gen - ref != span:
                l = f.readline()
                continue
            stat = toks[6].split('#')  # Pollak 0.01
            sampRes[-1].append(float(stat[1]))
            l = f.readline()
    ax.boxplot(sampRes)

plt.ioff()
for col in range(len(ncs)):
    for row in range(len(spans)):
        ax = plt.subplot(numRows, numCols, col * numRows + row + 1)
        doPlot(ax, ncs[col], spans[row], col == 0, row == numRows)
plt.savefig('fig1.png')
