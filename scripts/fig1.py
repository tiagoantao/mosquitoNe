import pylab
import sys
import ConfigParser

if len(sys.argv) not in [2]:
    print "syntax: %s <plotconfig>" % sys.argv[0]
    sys.exit(-1)

def getPlotConfig(fName):
    config = ConfigParser.ConfigParser()
    config.read(fName)
    return (eval(config.get("params","nc")), eval(config.get("params","span")),
            eval(config.get("params","sampleStrat")))


ncs, spans, sampleStrats = getPlotConfig(sys.argv[1])
numCols = len(ncs)
numRows = len(spans)

def doPlot(nc, span, startCol, endRow):
    sampRes = []
    for numIndivs, numLoci in sampleStrats:
        sampRes.append([])
        sampRes[(numIndivs, numLoci)] = []
        #we assume that t0 is file(gen[0])-1
        f = open("NeF-smp-%d-%d-%d.txt" % (numIndivs, numLoci, nc))
        #we assume span>1
        fileStartGen = int(f.readline().split(" ")[0])
        l = f.readline()
        while l != "":
            toks = l.rstrip().split(" ")
            gen = int(toks[0])
            if gen - fileStartGen + 1 == span:
                if toks[2]=="NeF":
                    sampRes[-1].append(float(toks[3]))
            l = f.readline()
    pylab.boxplot(sampRes)


for col in range(ncs):
    for row in range(spans):
        pylab.subplot(numRows, numCols, col*numRows + row + 1)
        doPlot(nc[col],spans[row], col==0, row==numRows)
