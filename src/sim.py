from simuPOP import *
from myUtils import getConfig
from sys import argv
import math
import numpy
import genepop
import os

if len(argv) not in [2, ]:
    print "Syntax:", argv[0], "<confFile>"
    exit(-1)
            
cfg = getConfig(argv[1])


def createGenome(size, numMSats, numSNPs):
    maxAlleleN = 100
    #print "Mutation model is most probably not correct", numMSats, numSNPs
    loci = (numMSats+numSNPs)*[1]
    initOps = []
    
    if numMSats > 0:
        diri = numpy.random.mtrand.dirichlet([1.0]*
                                             cfg.startAlleles, numMSats)
        if type(diri[0]) == float:
            diriList = diri
        else:
            diriList = list(diri[0])
                              
        initOps.append(
            InitGenotype(
                freq = [0.0] * ((maxAlleleN+1-8)//2) +
                    diriList +
                    [0.0] * ((maxAlleleN+1-8)//2),
                loci=range(0, numMSats)
           )
        )
    preOps = [StepwiseMutator(rates=cfg.mutFreq, loci=range(numMSats))]
    return loci, initOps, preOps
    

def createSinglePop(popSize, msatLoci, ins = None):
    preOps = [InitSex()]
    postOps = []
    pop = Population(size=[popSize], ploidy=2,
                     loci = [1]*msatLoci, chromTypes=[AUTOSOME]*msatLoci,)
    if ins:
        oExpr=('"samp/%s/%f/%%d/%%d/smp-%d-%d-%d-%%d-%%d.txt" %% ' +
            '(numIndivs, numLoci, gen, rep)')  % (ins, cfg.mutFreq, popSize, cfg.A, cfg.B)
    else:
        oExpr=('"samp/%f/%%d/%%d/smp-%d-%%d-%%d.txt" %% ' +
            '(numIndivs, numLoci, gen, rep)')  % (cfg.mutFreq, popSize)
    return pop, preOps, postOps, oExpr


def createSim(pop, reps):
    sim = Simulator (pop, rep = reps)
    return sim

def evolveSim(sim, gens, numMarkers, mateOp,
              genInitOps, genPreOps, popInitOps, popPostOps, reportOps,
              oExpr, saveGens):
    sim.evolve(
        initOps = genInitOps + popInitOps,
        preOps = genPreOps,
        postOps = popPostOps + reportOps,
        matingScheme = mateOp,
        gen=gens)
        


def saver(pop, param):
    oExpr = param
    for sampleStrat in cfg.sampleStrats:
        numIndivs, numLoci = sampleStrat
        rep = pop.dvars().rep
        gen = pop.dvars().gen
        ioExpr=eval(oExpr)
        genepop.save(pop, numIndivs, numLoci, ioExpr, 2)
    return True


def getSeasonFun(initialPop, genStart, A, B):
    def getSeasonFun(pop):
        gen = pop.dvars().gen 
        if gen<genStart:
            return initialPop
        else:
            return A*math.cos(2.0*(gen-genStart)*math.pi/12) + B
    return getSeasonFun

for sampleStrat in cfg.sampleStrats:
    numIndivs = sampleStrat[0]
    numLoci = sampleStrat[1]
    acum = os.getcwd()
    deepDirs = ['samp', cfg.mutFreq, numIndivs, numLoci]
    if cfg.demo == "bottle":
        deepDirs.insert(1, "bot")
    elif cfg.demo == "recover":
        deepDirs.insert(1, "rec")
    elif cfg.demo == "season":
        deepDirs.insert(1, "ses")
    for deepDir in deepDirs:
        try:
            os.mkdir(acum + os.sep + str(deepDir))
        except OSError:
            pass #OK
        acum += os.sep + str(deepDir)

if cfg.demo=="constant":
    (pop, popPreOps, popPostOps, oExpr) = createSinglePop(
        cfg.popSize, cfg.numMSats)
    (loci, genInitOps, genPreOps) = createGenome(
        cfg.popSize, cfg.numMSats, 0)
elif cfg.demo=="season":
    (pop, popPreOps, popPostOps, oExpr) = createSinglePop(
        cfg.popSize, cfg.numMSats, "ses")
    (loci, genInitOps, genPreOps) = createGenome(
        cfg.popSize, cfg.numMSats, 0)


if cfg.demo=="constant":
    mateOp = RandomMating()
elif cfg.demo=="season":
    mateOp = RandomMating(subPopSize=getSeasonFun(cfg.popSize,cfg.seasonGen,cfg.A,cfg.B))
    

reportOps = [
    PyOperator(func=saver, param = oExpr, at=cfg.saveGens),
    PyEval(r'"gen %d\n" % gen', reps=0),
    #PyEval(r'"size %s\n" % subPopSize', reps=0),
]
sim = createSim(pop, cfg.reps)
evolveSim(sim, cfg.gens, 0 + cfg.numMSats, mateOp,
              genInitOps, genPreOps, popPreOps, popPostOps, reportOps,
              oExpr, cfg.saveGens)
