import os
import sys
import shutil

from igrat.genetics.popgen.mlne.Controller import MLNEController as ctrl
import myUtils

if len(sys.argv) != 7:
    print("Syntax:", sys.argv[0], "<conffile> <rep> <ref> <gen> <inds> <loc>")
    sys.exit(-1)

etc = myUtils.getEtc()

cfg = myUtils.getConfig(sys.argv[1])
rep = int(sys.argv[2])
ref = int(sys.argv[3])
gen = int(sys.argv[4])
numIndivs = int(sys.argv[5])
numLoci = int(sys.argv[6])

mlnec = ctrl()


def calcMNE(rep, ref, gen, numIndivs, numLoci):
    gens = [ref, gen]
    fname = myUtils.getConc(cfg, numIndivs, numLoci, rep)
    tempName = "xaxa%d" % os.getpid()
    shutil.copyfile('%s-%d-%d' % (fname, ref, gen), tempName)
    point, ci = mlnec.run_mlne(open(tempName), False, True, 10000, gens)
    bottom, top = ci
    w = open('mne-%s-%d-%d-%d-%d-%d' % (sys.argv[1], rep, ref, gen, numIndivs, numLoci), 'w')
    w.write('%f %f %f\n' % (bottom, point, top))
    w.close()

calcMNE(rep, ref, gen, numIndivs, numLoci)
