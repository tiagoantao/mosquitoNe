from collections import defaultdict
import sys
#import ConfigParser

import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

import seaborn as sns

import myUtils
from myUtils import flt

ncs = [5000, 2000, 1000]
numIndivs = int(sys.argv[1])
numLoci = int(sys.argv[2])


def doLDNeCI(ax, nc, last_row):
    cfg = myUtils.getConfig('decl-%d-%d' % (nc, nc // 10))
    f = open(myUtils.getStatName(cfg, numIndivs, numLoci))
    l = f.readline()
    LD = [[], []]
    l = f.readline()
    gens = [cfg.declineGen, cfg.declineGen + 1, cfg.declineGen + 10]
    while l != 'temp\n':
        rep, gen = l.rstrip().split(' ')
        rep = int(rep)
        gen = int(gen)
        if gen not in gens:
            l = f.readline()
            l = f.readline()
            l = f.readline()
            l = f.readline()
            continue
        if gen == gens[0]:
            p = 0
        elif gen == gens[1]:
            p = 1
        else:
            p = 2
        l = f.readline()
        l = f.readline()
        toks = l.rstrip().split(' ')
        stat = toks[2].split('#')  # careful with pcrit
        if p != 2:  # We do not need the last result
            LD[p].append((flt(stat[0]), flt(stat[2])))
        l = f.readline()
        l = f.readline()
    bp_ld = []
    for i in range(len(LD)):
        p = 0 if i == 0 else 1
        bp_ld.append([y[p] if y[p] > 0 else 100000 for y in LD[i]])
    sns.boxplot(bp_ld, sym='', ax=ax)
    ax.set_xticklabels(['Bot CI\nNe=%d' % nc, 'Top CI\nNe=%d' % (nc / 10)])
    ax.set_ylim(0, nc // 5)


def doTempCI(ax, nc, last_row):
    cfg = myUtils.getConfig('decl-%d-%d' % (nc, nc // 10))
    f = open(myUtils.getStatName(cfg, numIndivs, numLoci))
    l = f.readline()
    l = f.readline()
    #study_gen = cfg.declineGen + 1
    study_ref = 49
    study_gens = [cfg.declineGen + 1, cfg.declineGen + 6,
                  cfg.declineGen + 12, cfg.declineGen + 24]
    while l != 'temp\n':
        l = f.readline()
    f.readline()  # pcrits
    l = f.readline()
    tdecl = defaultdict(list)
    while l != "":
        toks = l.rstrip().split(" ")
        #rep = int(toks[0])
        ref = int(toks[1])
        gen = int(toks[2])
        #if gen != study_gen:
        #    l = f.readline()
        #    continue
        if ref != study_ref:
            l = f.readline()
            continue
        if gen not in study_gens:
            l = f.readline()
            continue
        stat = toks[-1].split('#')  # Nei/Tajima 0+
        #tdecl[ref].append((flt(stat[0]), flt(stat[2])))
        tdecl[gen].append((flt(stat[0]), flt(stat[2])))
        l = f.readline()
    bp_temp = []
    refGens = sorted(list(tdecl.keys()))
    #for g in refGens:
    for g in study_gens:
        p = 1
        bp_temp.append([y[p] if y[p] > 0 and y[p] < 100000 else 100000 for y in tdecl[g]])
    if nc == 2000:
        print(bp_temp[1])
    sns.boxplot(bp_temp, sym='', ax=ax)
    ax.set_xticklabels([str(x - cfg.declineGen) for x in refGens])
    ax.set_ylim(0, nc // 1)


def doPlot(ax, nc, last_row):
    cfg = myUtils.getConfig('decl-%d-%d' % (nc, nc // 10))
    try:
        if nc == 1000:
            mdecl = defaultdict(dict)
            for rep, ref, gen, vals in myUtils.getMLNE(cfg, numIndivs, numLoci):
                top, point, bot = vals
                mdecl[ref].setdefault(gen, [])
                mdecl[ref][gen].append(point)
    except IOError:
        pass
    f = open(myUtils.getStatName(cfg, numIndivs, numLoci))
    l = f.readline()
    LD = defaultdict(list)
    coanc = defaultdict(list)
    het = defaultdict(list)
    l = f.readline()
    while l != 'temp\n':
        rep, gen = l.rstrip().split(' ')
        rep = int(rep)
        gen = int(gen)
        if gen < cfg.declineGen:
            l = f.readline()
            l = f.readline()
            l = f.readline()
            l = f.readline()
            continue
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
        if gen < cfg.declineGen:
            l = f.readline()
            continue
        tdecl[ref].setdefault(gen, [])
        stat = toks[-1].split('#')  # Nei/Tajima 0+
        tdecl[ref][gen].append(flt(stat[1]))
        l = f.readline()
    ax.plot([np.median([y if y > 0 else 100000 for y in LD[gen]])
             for gen in LD], label='LD')
    for ref, lst in tdecl.items():
        ax.plot([np.median([y if y > 0 else 100000 for y in lst[gen]])
                 for gen in lst], label=str(ref))
    if nc == 1000:
        for ref, lst in mdecl.items():
            ax.plot([np.median([y if y > 0 else 100000 for y in lst[gen]])
                    for gen in lst], label='m ' + str(ref))
    ax.set_xlim(0, 20)
    ax.set_ylim(0, nc)
    ax.legend()

plt.ioff()
numRows = len(ncs)
fig = plt.figure(figsize=(16, 9))
for row in range(numRows):
    ax = fig.add_subplot(numRows, 2, 1 + 2 * row)
    doPlot(ax, ncs[row], row == numRows - 1)
    ax = fig.add_subplot(numRows, 4, 3 + 4 * row)
    doLDNeCI(ax, ncs[row], row == numRows - 1)
    ax = fig.add_subplot(numRows, 4, 4 + 4 * row)
    doTempCI(ax, ncs[row], row == numRows - 1)
plt.savefig('decline-%d-%d.png' % (numIndivs, numLoci))
plt.savefig('decline-%d-%d.eps' % (numIndivs, numLoci))
