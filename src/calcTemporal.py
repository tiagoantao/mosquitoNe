import os
import sys
from copy import deepcopy
from scipy.stats import chi2
from myUtils import getConfig
from Bio.PopGen import GenePop

if len(sys.argv)!=2:
    print "Syntax:", sys.argv[0], "<confFile>"
    sys.exit(-1)

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
    size = len(list)
    if size % 2 == 0:
        return (list[size/2-1]+list[size/2])/2
    else:
        return list[size/2]

def getRefs(expr, gen, numReps):
    refs = []
    for rep in range(numReps):
        fname = eval(expr)
        fi = open(fname)
        rec = GenePop.parse(fi)
        fi.close()
        csld = GenePop.Record()
        csld.loci_list = deepcopy(rec.loci_list)
        csld.populations=[[]]
        for pop in rec.populations:
            csld.populations[0].extend(pop)
        refs.append(csld)
    return refs

def calcNe(F, genDiff, pop1Len, pop2Len):
    fk, df = F
    Nes = []
    for ft in [fk, df*fk/chi2.isf(0.025,df), df*fk/chi2.isf(0.975, df)]:
        Ne = (1.0*(genDiff))/ (2*(
             ft
               - 1.0/(2*pop1Len)
               - 1.0/(2*pop2Len)
             ))
        if Ne<0 or Ne>10000: Ne = float('inf')
        Nes.append(Ne)
    return Nes

def calcFs(expr, numReps, gpRefs, gen, refGen):
    for rep in range(numReps):
        fname = eval(expr)
        gpRef = gpRefs[rep]
        fi = open(fname)
        rec = GenePop.parse(fi)
        fi.close()
        Nes = []
        NesP = []
        for pop in rec.populations:
            F = calcF(gpRef.populations[0], pop, len(gpRef.loci_list))
            Mes = calcNe(F, gen-refGen, len(gpRef.populations[0]), len(pop))
            Ne = Mes[0]
            NeP = Mes[1:]
            #print Ne
            #print "F", F, 100.0*F/129.0, 100*F/74.0
            #Ne = (1.0*(gen-refGen))/ (2*F)
            NesP.append(str(min(NeP))+"#"+str(max(NeP)))
            Nes.append(Ne)
        print gen, rep, "NeF", " ".join(map(lambda x: str(x), Nes))
        print gen, rep, "NeFPow", " ".join(map(lambda x: str(x), NesP))


def calcF(refPop, pop, numLoci):
    Fk = []
    for i in range(numLoci):
        refList = {}
        popList = {}
        cnt = 0
        for indiv in refPop:
            popList.setdefault(indiv[1][i][0],0)
            popList.setdefault(indiv[1][i][1],0)
            refList[indiv[1][i][0]] = refList.get(indiv[1][i][0],0)+1
            refList[indiv[1][i][1]] = refList.get(indiv[1][i][1],0)+1
            cnt += 2
        cnt = 0
        for indiv in pop:
            refList.setdefault(indiv[1][i][0],0)
            refList.setdefault(indiv[1][i][1],0)
            popList[indiv[1][i][0]] = popList.get(indiv[1][i][0],0)+1
            popList[indiv[1][i][1]] = popList.get(indiv[1][i][1],0)+1
            cnt += 2
        for key in refList.keys():
            refList[key] = 1.0*refList[key]/cnt
            popList[key] = 1.0*popList[key]/cnt
        sum = 0.0
        for all in refList.keys():
            if refList[all]>0.0 or popList[all]>0.0:
                sum += ((refList[all] - popList[all])**2)/(
                        ((refList[all] + popList[all])/2.0)
                       )
        numAll = len(refList.keys())
        if numAll>1:
            Fk.append((numAll,(1.0/(numAll-1.0)*sum)))
    numerator = 0.0
    denominator = 0.0
    for na, fk in Fk:
        denominator += na - 1
        numerator += (na - 1) * fk
    F = numerator / denominator
    # denominator is the number of degrees of freedom
    return F, denominator


for refGen in cfg.refGens:
    for sampleStrat in cfg.sampleStrats:
        numIndivs, numLoci = sampleStrat
        if cfg.demo=="constant":
            oExpr=('"samp/%f/%d/%d/smp-%d-%%d-%%d.txt" %% (gen, rep,)'  %
               (cfg.mutFreq, numIndivs, numLoci, cfg.popSize) )
        elif cfg.demo=="season":
            oExpr=('"samp/ses/%f/%d/%d/smp-%d-%d-%d-%%d-%%d.txt" %% (gen, rep,)'  %
               (cfg.mutFreq, numIndivs, numLoci, cfg.popSize, cfg.A, cfg.B) )

        gpRefs = getRefs(oExpr, refGen, cfg.reps)

        old = sys.stdout
        if cfg.demo=="constant":
            sys.stdout = open("NeF-smp-%d-%d-%d-%d.txt" %
                          (numIndivs, numLoci, cfg.popSize, refGen), "w")
        elif cfg.demo=="season":
            sys.stdout = open("NeF-ses-%d-%d-%d-%d-%d-%d.txt" %
                          (numIndivs, numLoci, cfg.popSize, refGen, cfg.A, cfg.B), "w")

        for gen in cfg.futureGens:
            calcFs(oExpr, cfg.reps, gpRefs, gen, refGen)
        stdout = old
