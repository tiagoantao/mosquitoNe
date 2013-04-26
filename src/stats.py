import os
import sys
import shutil
from copy import deepcopy

from Bio.PopGen.GenePop.Controller import GenePopController
from PopGen import NeEstimator2
from PopGen.NeEstimator2.Controller import NeEstimator2Controller
import myUtils

if len(sys.argv) != 2:
    print "Syntax:", sys.argv[0], "<conffile>"
    sys.exit(-1)

etc = myUtils.getEtc()

gpc = GenePopController(etc['genepop'])
ne2 = NeEstimator2Controller(etc['ne2'])

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
    blocks = myUtils.getBlocks(cfg)
    for rep in range(cfg.reps):
        fname = myUtils.getConc(cfg, numIndivs, numLoci, rep)
        tempName = "xaxa%d" % os.getpid()
        coanc = {}
        ldne = {}
        hetNb = {}
        temp = []
        for i in range(len(blocks)):
            gens = blocks[i]
            shutil.copyfile(fname + "-" + str(i), tempName)
            ne2.run_neestimator2(tempName, tempName + '.out',
                                 LD=True, hets=True, coanc=True,
                                 temp=gens)
            ldout = open(tempName + '.out')
            rec = NeEstimator2.parse(ldout)
            for j in range(len(gens)):
                if gens[j] in cfg.futureGens:
                    coanc[gens[j]] = rec.coanc[j]["EstNeb"]
                    ldne[gens[j]] = [rec.ld[j][c]["EstNe"] for c in
                                     range(len(rec.ld[j]))]
                    hetNb[gens[j]] = [rec.het[j][c]["EstNeb"] for c in
                                      range(len(rec.het[j]))]
            for tmp in rec.temporal:
                temp.append(tmp)
        for gen in cfg.futureGens:
            print rep, gen,
            print coanc[gen]
            for ld in ldne[gen]:
                print ld,
            print
            for hnb in hetNb[gen]:
                print hnb,
            print
        for tmp in temp:
            g1 = tmp["generation1"]
            g2 = tmp["generation2"]
            if g1 in cfg.refGens and g2 in cfg.futureGens:
                pl = tmp["results"]["Pollak"]
                jr = tmp["results"]["Jorde/Ryman"]
                nt = tmp["results"]["Nei/Tajima"]
                print rep, "temp", g1, g2,
                print "#".join([str(pl[c]["Ne"]) for c in range(len(pl))]),
                print "#".join([str(jr[c]["Ne"]) for c in range(len(pl))]),
                print "#".join([str(nt[c]["Ne"]) for c in range(len(pl))])

stdout = sys.stdout
for numIndivs, numLoci in cfg.sampleStrats:
    sys.stdout = open(myUtils.getStatName(cfg, numIndivs, numLoci), "w")
    calcStats(numIndivs, numLoci)
sys.stdout = stdout
