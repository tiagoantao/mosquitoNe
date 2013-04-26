import math
import os
import sys

import pylab

import myUtils

if len(sys.argv) != 2:
    print "Syntax:", sys.argv[0], "<conffile>"
    sys.exit(-1)

etc = myUtils.getEtc()

cfg = myUtils.getConfig(sys.argv[1])

for t in range(cfg.seasonGen, cfg.gens):
    y = cfg.A * math.cos(2 * math.pi * (t - cfg.seasonGen) / cfg.T) + cfg.B
    print t, y
sys.exit(0)

for numIndivs, numLoci in cfg.sampleStrats:
    fname = myUtils.getStatName(cfg, numIndivs, numLoci)
    for rec in myUtils.getStat(open(fname)):
        print rec
