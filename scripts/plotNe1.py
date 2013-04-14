import math
import pylab
from myUtils import getConfig
import sys

if len(sys.argv) not in [2,3]:
    print "syntax: %s <config> [top|bottom]" % sys.argv[0]
    sys.exit(-1)

cfg = getConfig(sys.argv[1])
myFun=float
myStat = "NeF"
if len(sys.argv) == 3:
    myStat = "NeFPow"
    if sys.argv[2] == "top":
        myFun = lambda x: float(x.split("#")[1])

def getStat(name, f, fun=float):
    l = f.readline()
    while l!="":
        toks = l.rstrip().split(" ")
        if toks[2]==name:
            return int(toks[0]), int(toks[1]), fun(toks[3])
        l = f.readline()
    return None

def getNc(gen):
    #This is cut n paste (sim.py)
    return cfg.A*math.cos(2.0*(gen-cfg.seasonGen)*math.pi/12) + cfg.B

for sampleStrat in cfg.sampleStrats:
    for refGen in cfg.refGens:
        pylab.clf()
        numIndivs = sampleStrat[0]
        numLoci = sampleStrat[1]
        vals = {}
        if cfg.demo=="constant":
            f = open("NeF-smp-%d-%d-%d.txt" % (numIndivs, numLoci, cfg.popSize))
        elif cfg.demo=="season":
            f = open("NeF-ses-%d-%d-%d-%d-%d-%d.txt" % (numIndivs, numLoci, cfg.popSize, refGen, cfg.A, cfg.B))
        try:
            while True:
                gen, iter, nef = getStat(myStat, f, myFun)
                pos = gen - refGen
                lst = vals.setdefault(pos,[])
                lst.append(nef)
        except TypeError:
            pass #OK
        f.close()
        keys = vals.keys()
        keys.sort()
        x = []
        y975 = []
        y50 = []
        y025 = []
        nc = []
        for key in keys:
            x.append(key)
            nc.append(getNc(key))
            y = []
            for val in vals[key]:
                    y.append(val)
            y.sort()
            try:
                y975.append(y[975*len(y)/1000])
            except ZeroDivisionError:
                y975.append(None)
            try:
                y50.append(float(len(y))/reduce(lambda x,y:1.0/y+x,y,0))
            except TypeError:
                y50.append(0.0)
            except ZeroDivisionError:
                y50.append(None)
            try:
                y025.append(y[25*len(y)/1000])
            except ZeroDivisionError:
                y025.append(None)

        pylab.plot(x,y975,label="90")
        pylab.plot(x,y50,label="50")
        pylab.plot(x,y025,label="10")
        pylab.plot(x,nc,label="Nc")
        if cfg.demo=="constant":
            pylab.title("size: %d - inds: %d - msats: %d" % (cfg.popSize, numIndivs, numLoci))
        if cfg.demo=="season":
            pylab.title("size: %d - inds: %d - msats: %d - refGen: %d" % (cfg.popSize, numIndivs, numLoci, refGen))
        v=list(pylab.axis())
        v[0]=0
        v[1]=40
        v[3]=2400
        pylab.axis(v)
        #pylab.show()
        pylab.legend()
        pylab.savefig("%s-%d-%d-%d-%d.eps" % (myStat, cfg.popSize, numIndivs, numLoci, refGen))

