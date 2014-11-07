import ConfigParser
import sys

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

import seaborn as sns

import myUtils
from myUtils import flt

if len(sys.argv) != 2:
    print "Syntax:", sys.argv[0], "<conffile>"
    sys.exit(-1)

etc = myUtils.getEtc()


def getPlotConfig(fName):
    config = ConfigParser.ConfigParser()
    config.read(fName)
    return (eval(config.get("params", "nc")),
            eval(config.get("params", "span")),
            eval(config.get("params", "sampleStrat")))

ncs, spans, sampleStrats = getPlotConfig(sys.argv[1])


def plotComparison(ax, nc, sampleStrat, startCol, endRow):
    cfg = myUtils.getConfig('simple%d' % nc)
    numIndivs, numLoci = sampleStrat
    if numIndivs == 60 and numLoci == 20 and nc == 1000:
        do_mlne = True
    else:
        do_mlne = False
    if do_mlne:
        myml = {}
        for span in spans:
            myml[span] = []
        for rep, ref, gen, vals in myUtils.getMLNE(cfg, numIndivs, numLoci):
            dist = gen - ref
            if dist not in spans:
                continue
            val = vals[1]
            myml[dist].append(val)
    fname = myUtils.getStatName(cfg, numIndivs, numLoci)
    f = open(fname)
    f.readline()
    l = f.readline()
    LD = {}
    for gen in cfg.futureGens:
        LD[gen] = []
    while l != 'temp\n':
        rep, gen = l.rstrip().split(' ')
        rep = int(rep)
        gen = int(gen)
        if gen not in cfg.futureGens:
            continue
        l = f.readline()
        #toks = l.rstrip().split(' ')
        #coanc[gen].append(flt(toks[1]))
        l = f.readline()
        toks = l.rstrip().split(' ')
        stat = toks[3].split('#')  # careful with pcrit
        st = flt(stat[1])
        LD[gen].append(st if st < 100000 else 100000)
        l = f.readline()
        #toks = l.rstrip().split(' ')
        #het[gen].append(flt(toks[1]))  # pcrit...
        l = f.readline()
    f.readline()  # pcrit spec
    l = f.readline()
    mytemp = {}
    for span in spans:
        mytemp[span] = []
    while l != '':
        toks = l.rstrip().split(' ')
        ref = int(toks[1])
        gen = int(toks[2])
        dist = gen - ref
        if dist not in spans:
            l = f.readline()
            continue
        stat = toks[-3].split('#')  # careful with pcrit
        st = flt(stat[1])
        mytemp[dist].append(st if st > 0 and st < 100000 else 100000)
        l = f.readline()
    sns.boxplot([mytemp[span] for span in spans] + [[LD[cfg.futureGens[0]]]],
                sym='', ax=ax)
    if do_mlne:
        sns.boxplot([myml[span] for span in spans] + [],
                    sym='', ax=ax, widths=.5, color='bright')
    ax.set_ylim(0, 3 * nc)
    #ax.get_yaxis().set_ticks([nc // 2, nc, 3 * nc // 2, 2 * nc])
    ax.get_yaxis().set_ticks([nc // 2, nc, 3 * nc // 2, 2 * nc, 3 * nc])
    ax.axhline(nc)
    if not startCol:
        ax.set_yticklabels(['', '', '', ''])
    if endRow:
        ax.set_xticklabels(['span-%d' % span for span in spans] + ['LD'])

plt.ioff()
numCols = len(sampleStrats)
numRows = len(ncs)
fig, axs = plt.subplots(numRows, numCols, sharex=True, figsize=(16, 9))
for row in range(len(ncs)):
    for col in range(len(sampleStrats)):
        ax = axs[row, col]
        plotComparison(ax, ncs[row], sampleStrats[col], col == 0, row ==
                       numRows - 1)
        if row == 0:
            ymin, ymax = ax.get_ylim()
            xmin, xmax = ax.get_xlim()
            ni, nl = sampleStrats[col]
            ax.text((xmax - xmin) / 2, ymax, 'I = %d, L = %d' % (ni, nl),
                    va='bottom', ha='center')
fig.tight_layout(h_pad=0, w_pad=0.1)
fig.savefig('comparison.png')
fig.savefig('comparison.eps')
