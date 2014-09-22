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

# for non-temporal
for numIndivs, numLoci in cfg.sampleStrats:
    for rep in range(cfg.reps):
        wname = myUtils.getConc(cfg, numIndivs, numLoci, rep)
        w = open(wname, "w")
        w.write("conc\n")
        gens = cfg.futureGens
        gens.sort()
        onFirst = True
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

# for temporal
for numIndivs, numLoci in cfg.sampleStrats:
    for rep in range(cfg.reps):
        for ref in cfg.refGens:
            rname = myUtils.getExpr(cfg, numIndivs, numLoci, ref, rep, True)
            f = open(rname)
            begin = f.readlines()
            for me in cfg.futureGens:
                wname = myUtils.getConc(cfg, numIndivs, numLoci, rep)
                w = open('%s-%d-%d' % (wname, ref, me), "w")
                w.write("conc\n")
                #w.write(''.join(begin[1:numLoci + 1]))
                w.write(''.join(begin[1: - numIndivs // 2]))
                name = myUtils.getExpr(cfg, numIndivs, numLoci, me, rep, True)
                f = open(name)
                onPop = False
                i = 0
                for l in f:
                    if onPop:
                        #w.write(begin[numLoci + i + 2])
                        w.write('t' + l)
                        i += 1
                        if i == numIndivs // 2:
                            break
                    elif l.startswith("Pop"):
                        w.write(l)
                        onPop = True
            w.close()
