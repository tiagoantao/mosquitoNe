import sys
from copy import deepcopy

from genomics.parallel.executor import Local

from igrat.genetics.popgen.mlne.Controller import MLNEController as ctrl
import myUtils

if len(sys.argv) != 2:
    print("Syntax:", sys.argv[0], "<conffile>")
    sys.exit(-1)

etc = myUtils.getEtc()

ex = Local(-40)
cfg = myUtils.getConfig(sys.argv[1])


mlnec = ctrl()


def calcMNE(numIndivs, numLoci):
    for rep in range(cfg.reps):
        for ref in cfg.refGens:
            for gen in cfg.futureGens:
                print('python3', 'src/single_mlne.py %s %d %d %d %d %d' % (
                    sys.argv[1], rep, ref, gen, numIndivs, numLoci))
                ex.submit('python3', 'src/single_mlne.py %s %d %d %d %d %d' % (
                    sys.argv[1], rep, ref, gen, numIndivs, numLoci))

stdout = sys.stdout
cfg.sampleStrats = [(60, 20)]
for numIndivs, numLoci in cfg.sampleStrats:
    #sys.stdout = open(myUtils.getStatName(cfg, numIndivs, numLoci) + '-mne', "w")
    calcMNE(numIndivs, numLoci)
sys.stdout = stdout
