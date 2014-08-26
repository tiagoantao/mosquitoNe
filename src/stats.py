import os
import sys
import shutil
from copy import deepcopy

from Bio.PopGen.GenePop.Controller import GenePopController
from igrat.genetics.popgen import ne2
from igrat.genetics.popgen.ne2.Controller import NeEstimator2Controller
import myUtils
from myUtils import flt

if len(sys.argv) != 2:
    print "Syntax:", sys.argv[0], "<conffile>"
    sys.exit(-1)

etc = myUtils.getEtc()

gpc = GenePopController(etc['genepop'])
ne2c = NeEstimator2Controller(etc['ne2'])

cfg = myUtils.getConfig(sys.argv[1])


def all(list):
    for e in list:
        print e,
    print


def interval(list):
    for e1, e2 in list:
        print str(e1) + "#" + str(e2),
    print


def median(list):
    list = deepcopy(list)
    list.sort()
    size = len(list)
    if size % 2 == 0:
        return (list[size / 2 - 1] + list[size / 2]) / 2
    else:
        return list[size / 2]


def calcStats(numIndivs, numLoci):
    crits = [0.05, 0.02, 0.01, 0]
    print ' '.join([str(x) for x in crits])
    for rep in range(cfg.reps):
        fname = myUtils.getConc(cfg, numIndivs, numLoci, rep)
        tempName = "xaxa%d" % os.getpid()
        coanc = {}
        ldne = {}
        ldneCI = {}
        hetNb = {}
        shutil.copyfile(fname, tempName)
        ne2c.run_neestimator2('.', tempName, '.', tempName + '.out',
                              crits=crits, LD=True, hets=True, coanc=True)
        ldout = open(tempName + '.out')
        rec = ne2.parse(ldout)
        for j, gen in enumerate(cfg.futureGens):
            coanc[gen] = rec.coanc[j]["EstNeb"]
            ldne[gen] = [rec.ld[j][c]["EstNe"] for c in range(len(rec.ld[j]))]
            ldneCI[gen] = [rec.ld[j][c]["ParaNe"] for c in range(len(rec.ld[j]))]
            hetNb[gen] = [rec.het[j][c]["EstNeb"] for c in range(len(rec.het[j]))]
            print rep, gen
            print 'coanc',
            print coanc[gen]
            print 'ld',
            for i, ld in enumerate(ldne[gen]):
                low, high = ldneCI[gen][i]
                print '%f#%f#%f' % (flt(low), flt(ld), flt(high)),
            print
            print 'het',
            for hnb in hetNb[gen]:
                print hnb,
            print


def calcTemporalStats(numIndivs, numLoci):
    crits = [0.05, 0.02, 0.01, 0]
    print 'temp'
    print ' '.join([str(x) for x in crits])
    for rep in range(cfg.reps):
        for ref in cfg.refGens:
            for gen in cfg.futureGens:
                gens = [(ref, gen)]
                fname = myUtils.getConc(cfg, numIndivs, numLoci, rep)
                tempName = "xaxa%d" % os.getpid()
                for me in cfg.futureGens:
                    shutil.copyfile('%s-%d-%d' % (fname, ref, me), tempName)
                    ne2c.run_neestimator2('.', tempName, '.', tempName + '.out',
                                          LD=False, crits=crits, temp=gens)
                    ldout = open(tempName + '.out')
                    rec = ne2.parse(ldout)
                    print rep, ref, me,
                    case = rec.temporal[0]
                    for method, results in case['results'].items():
                        print method,
                        for crit in results:
                            print '%.1f#%.1f#%.1f' % (
                                flt(crit['ParaTemp'][0]), flt(crit['Ne']),
                                flt(crit['ParaTemp'][1])),
                    print

stdout = sys.stdout
for numIndivs, numLoci in cfg.sampleStrats:
    sys.stdout = open(myUtils.getStatName(cfg, numIndivs, numLoci), "w")
    calcStats(numIndivs, numLoci)
    calcTemporalStats(numIndivs, numLoci)
sys.stdout = stdout
