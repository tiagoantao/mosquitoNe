import os
import sys

import myUtils

model = sys.argv[1]

try:
    os.mkdir('samp/conc')
except OSError:
    pass  # Alread exists: OK


cfg = myUtils.getConfig(model)
popSize = cfg.popSize
blocks = myUtils.getBlocks(cfg)

for numIndivs, numLoci in cfg.sampleStrats:
    for rep in range(cfg.reps):
        for i in range(len(blocks)):
            wname = myUtils.getConc(cfg, numIndivs, numLoci, rep)
            w = open(wname + "-" + str(i), "w")
            onFirst = True
            w.write("conc\n")
            gens = set(blocks[i])
            gens = list(gens)
            gens.sort()
            for gen in gens:
                name = myUtils.getExpr(cfg, numIndivs, numLoci, gen, rep, True)
                f = open(name)
                sendOut = False
                f.readline()
                for l in f:
                    if sendOut or onFirst:
                        w.write(l)
                    elif l.startswith("Pop"):
                        w.write(l)
                        sendOut = True
                onFirst = False
            w.close()
