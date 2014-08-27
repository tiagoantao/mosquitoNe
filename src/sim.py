import sys
import math
import numpy
import os

import simuPOP as sp
import genepop
import myUtils

if len(sys.argv) not in [2, ]:
    print "Syntax:", sys.argv[0], "<confFile>"
    exit(-1)

cfg = myUtils.getConfig(sys.argv[1])


def createGenome(size, numMSats, numSNPs):
    maxAlleleN = 100
    #print "Mutation model is most probably not correct", numMSats, numSNPs
    loci = (numMSats + numSNPs) * [1]
    initOps = []

    if numMSats > 0:
        diri = numpy.random.mtrand.dirichlet([1.0] *
                                             cfg.startAlleles, numMSats)
        if type(diri[0]) == float:
            diriList = diri
        else:
            diriList = list(diri[0])

        initOps.append(
            sp.InitGenotype(
                freq=[0.0] * ((maxAlleleN + 1 - 8) // 2) +
                diriList + [0.0] * ((maxAlleleN + 1 - 8) // 2),
                loci=range(0, numMSats))
        )
    if cfg.mutFreq > 0:
        preOps = [sp.StepwiseMutator(rates=cfg.mutFreq, loci=range(numMSats))]
    else:
        preOps = []
    return loci, initOps, preOps


def createSinglePop(popSize, msatLoci, model):
    preOps = [sp.InitSex()]
    postOps = []
    pop = sp.Population(size=[popSize], ploidy=2, loci=[1] * msatLoci,
                        chromTypes=[sp.AUTOSOME] * msatLoci,)
    oExpr = myUtils.getExpr(cfg, numIndivs, numLoci, None, None, False)
    return pop, preOps, postOps, oExpr


def createSim(pop, reps):
    sim = sp.Simulator(pop, rep=reps)
    return sim


def evolveSim(sim, gens, numMarkers, mateOp,
              genInitOps, genPreOps, popInitOps, popPostOps, reportOps,
              oExpr, saveGens):
    sim.evolve(
        initOps=genInitOps + popInitOps,
        preOps=genPreOps,
        postOps=popPostOps + reportOps,
        matingScheme=mateOp,
        gen=gens)


def saver(pop, param):
    oExpr = param
    for sampleStrat in cfg.sampleStrats:
        numIndivs, numLoci = sampleStrat
        rep = pop.dvars().rep
        gen = pop.dvars().gen
        ioExpr = eval(oExpr)
        genepop.save(pop, min([numIndivs, pop.popSize()]), numLoci, ioExpr, 2)
    return True


def getSeasonFun(initialPop, genStart, A, B, T):
    def getSeasonFun(pop):
        gen = pop.dvars().gen
        if gen < genStart:
            return initialPop
        else:
            return A * math.cos(2.0 * (gen - genStart) * math.pi / T) + B
    return getSeasonFun


def getDeclineFun(initialPop, declineGen, smallPop):
    def getDeclineFun(pop):
        gen = pop.dvars().gen
        if gen < declineGen:
            return initialPop
        else:
            return smallPop
    return getDeclineFun


for sampleStrat in cfg.sampleStrats:
    numIndivs = sampleStrat[0]
    numLoci = sampleStrat[1]
    acum = os.getcwd()
    deepDirs = ['samp', "%f" % cfg.mutFreq, numIndivs, numLoci]
    if cfg.demo == "decline":
        deepDirs.insert(1, "decl")
    elif cfg.demo == "recover":
        deepDirs.insert(1, "rec")
    elif cfg.demo == "season":
        deepDirs.insert(1, "ses")
    for deepDir in deepDirs:
        try:
            os.mkdir(acum + os.sep + str(deepDir))
        except OSError:
            pass  # OK
        acum += os.sep + str(deepDir)

pop, popPreOps, popPostOps, oExpr = createSinglePop(
    cfg.popSize, cfg.numMSats, cfg.demo)
loci, genInitOps, genPreOps = createGenome(cfg.popSize, cfg.numMSats, 0)


if cfg.demo == "constant":
    mateOp = sp.RandomMating()
elif cfg.demo == "season":
    mateOp = sp.RandomMating(subPopSize=getSeasonFun(cfg.popSize,
                                                     cfg.seasonGen,
                                                     cfg.A, cfg.B, cfg.T))
elif cfg.demo == "decline":
    mateOp = sp.RandomMating(subPopSize=getDeclineFun(cfg.popSize,
                                                      cfg.declineGen,
                                                      cfg.declineSize))

reportOps = [
    sp.PyOperator(func=saver, param=oExpr, at=cfg.saveGens),
    sp.Stat(popSize=True),
    sp.PyEval(r'"gen %d" % gen', reps=0),
    sp.PyEval(r'" size %d\n" % subPopSize[0]', reps=0),
]
sim = createSim(pop, cfg.reps)
evolveSim(sim, cfg.gens, 0 + cfg.numMSats, mateOp,
          genInitOps, genPreOps, popPreOps, popPostOps, reportOps,
          oExpr, cfg.saveGens)
