from Bio.PopGen.GenePop.Controller import GenePopController
from Bio.PopGen.LDNe.Controller import LDNeController
from Bio.PopGen import LDNe
import os
import sys
import shutil
from copy import deepcopy
from myUtils import getConfig

if len(sys.argv)!=2:
    print "Syntax:", sys.argv[0], "<conffile>"
    sys.exit(-1)


gpc = GenePopController('/home/tiago/bio')
ldnec = LDNeController('/home/tiago/bio')

cfg = getConfig(sys.argv[1])


def all(list):
    for e in list:
        print e,
    print

def interval(list):
    for e1,e2 in list:
        print str(e1) + "#" + str(e2),
    print

def median(list):
    list = deepcopy(list)
    list.sort()
    size = len(list)
    if size % 2 == 0:
        return (list[size/2-1]+list[size/2])/2
    else:
        return list[size/2]

def calcStats(numIndivs, numLoci):
    for gen in cfg.saveGens:
        for rep in range(cfg.reps):
           #print rep,
           if cfg.demo=="constant":
                fname=('samp/%f/%d/%d/smp-%d-%d-%d.txt'  %
                   (cfg.mutFreq, numIndivs, numLoci, cfg.popSize, gen, rep) )
           tempName =  "xaxa%d" % (os.getpid(),)
           shutil.copyfile(fname, tempName)
           ldnec.run_ldne(tempName, tempName + '.out')
           ldout = open(tempName + '.out')
           ldres = LDNe.RecordParser().parse(ldout)
           mNes = []
           mNesPow = []
           for id, fcases in ldres.populations:
               if numIndivs<25:
                   hm, ic, or2, er2, (ne, (ne95, ne05), (j95, j05)) = fcases[0]
               else:
                   hm, ic, or2, er2, (ne, (ne95, ne05), (j95, j05)) = fcases[1]
               mNes.append(ne)
               mNesPow.append((ne95,ne05))
           ldout.close()

           genos =  gpc.get_loci_genotype_counts(tempName)
           mGenos = []
           mRich = []
           for popGenos in genos[0]:
               currGenos = []
               currRich = []
               for loci, lst, expHo, obsHo, expHe, obsHe in popGenos:
                   alleles = set()
                   for a1, a2, c in lst:
                       alleles.add(a1)
                       alleles.add(a2)
                   currGenos.append(expHe/(expHe+expHo))
                   currRich.append(len(alleles))
               mGenos.append(currGenos)
               mRich.append(currRich)

           for pi in range(len(mGenos)):
               print gen, rep, "ExpHe" + str(pi+1),
               all(mGenos[pi])
               print gen, rep, "AllRich" + str(pi+1),
               all(mRich[pi])
           print gen, rep, "Ne",
           all(mNes)
           print gen, rep, "NePow",
           interval(mNesPow)

stdout = sys.stdout
for numIndivs, numLoci in cfg.sampleStrats:
    if cfg.demo=="constant":
        sys.stdout = open("LDNe-smp-%d-%d-%d.txt" %
                          (numIndivs, numLoci, cfg.popSize), "w")
    calcStats(numIndivs, numLoci)
sys.stdout=out
