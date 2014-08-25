import ConfigParser
import os

sysDir = "/".join(os.path.realpath(__file__).split("/")[:-2])
confDir = sysDir + "/conf/"
etcConf = sysDir + "/etc.conf"

oExpr = {
    "constant": '"samp/%f/%%d/%%d/smp-%d-%%d-%%d.txt" %% ',
    "season":   '"samp/ses/%f/%%d/%%d/smp-%d-%d-%d-%d-%%d-%%d.txt" %% ',
    "bot":      '"samp/bot/%f/%%d/%%d/smp-%d-%d-%d-%%d-%%d.txt" %% ',
}


def flt(x):
    try:
        return float(x)
    except ValueError:
        return float("inf")


def getEtc():
    etc = {}
    cfg = ConfigParser.ConfigParser()
    cfg.read(etcConf)
    etc["ne2"] = cfg.get("bin", "ne2")
    etc["genepop"] = cfg.get("bin", "genepop")
    return etc


def getConfig(fName):
    config = ConfigParser.ConfigParser()
    config.read(confDir + fName + ".conf")

    class Cfg:
        pass

    cfg = Cfg()
    cfg.demo = config.get("pop", "demo")
    if cfg.demo == "constant":
        cfg.popSize = config.getint("pop", "popSize")
    if cfg.demo == "season":
        cfg.popSize = config.getint("pop", "popSize")
        cfg.seasonGen = config.getint("pop", "seasonGen")
        cfg.A = config.getint("pop", "A")
        cfg.B = config.getint("pop", "B")
        try:
            cfg.T = config.getint("pop", "T")
        except ConfigParser.NoOptionError:
            cfg.T = 12

    cfg.startAlleles = config.getint("genome", "startAlleles")
    cfg.mutFreq = config.getfloat("genome", "mutFreq")
    cfg.numMSats = config.getint("genome", "numMSats")

    cfg.reps = config.getint("sim", "reps")
    cfg.gens = config.getint("sim", "gens")
    cfg.sampleStrats = eval(config.get("sim", "sampleStrats"))
    cfg.saveGens = eval(config.get("sim", "saveGens"))

    cfg.refGens = eval(config.get("stat", "refGens"))
    cfg.futureGens = eval(config.get("stat", "futureGens"))
    return cfg


def getExpr(cfg, numIndivs, numLoci, gen, rep, doeval=False):
    if cfg.demo == "season":
        rExpr = (oExpr[cfg.demo] +
                 '(numIndivs, numLoci, gen, rep)') % (
                     cfg.mutFreq, cfg.popSize, cfg.A, cfg.B, cfg.T)
    elif cfg.demo == "bottle":
        rExpr = (oExpr[cfg.demo] +
                 '(numIndivs, numLoci, gen, rep)') % (
                     cfg.mutFreq, cfg.popSize, cfg.bottleGen, cfg.popSize2)
    elif cfg.demo == "constant":
        rExpr = (oExpr[cfg.demo] +
                 '(numIndivs, numLoci, gen, rep)') % (cfg.mutFreq, cfg.popSize)
    if doeval:
        return eval(rExpr)
    else:
        return rExpr


def getConc(cfg, numIndivs, numLoci, rep):
    if cfg.demo == "season":
        return "samp/conc/ses-%d-%d-%f-%d-%d-%d-%d-%d.txt" % (
            numIndivs, numLoci, cfg.mutFreq, cfg.popSize,
            cfg.A, cfg.B, cfg.T, rep)
    elif cfg.demo == "bottle":
        return "samp/conc/bot-%d-%d-%f-%d-%d-%d-%d.txt" % (
            numIndivs, numLoci, cfg.mutFreq, cfg.popSize,
            cfg.bottleGen, cfg.popSize2, rep)
    elif cfg.demo == "constant":
        return "samp/conc/con-%d-%d-%f-%d-%d.txt" % (
            numIndivs, numLoci, cfg.mutFreq, cfg.popSize, rep)


def getStatName(cfg, numIndivs, numLoci):
    if cfg.demo == "constant":
        return "con-%d-%d-%d.txt" % (numIndivs, numLoci, cfg.popSize)
    elif cfg.demo == "season":
        return "ses-%d-%d-%d-%d-%d-%d.txt" % (numIndivs, numLoci,
                                              cfg.popSize, cfg.A,
                                              cfg.B, cfg.T)
    elif cfg.demo == "bottle":
        return "ses-%d-%d-%d-%d-%d.txt" % (numIndivs, numLoci,
                                           cfg.popSize, cfg.bottleGen,
                                           cfg.popSize2)


def getStat(f):
    l = f.readline()
    while l != "":
        rec = {}
        l = l.rstrip()
        if l.find("temp") > -1:
            toks = l.split(" ")
            rec["type"] = "temp"
            rec["rep"] = int(toks[0])
            rec["g1"] = int(toks[2])
            rec["g2"] = int(toks[3])
            rec["Pollak"] = [flt(x) for x in toks[4].split("#")]
            rec["Jorde/Ryman"] = [flt(x) for x in toks[5].split("#")]
            rec["Nei/Tajima"] = [flt(x) for x in toks[6].split("#")]
            yield rec
        else:
            toks = l.split(" ")
            rec["type"] = "notemp"
            rec["rep"] = int(toks[0])
            rec["gen"] = int(toks[1])
            rec["coanc"] = flt(toks[2])
            rec["LD"] = [flt(x) for x in f.readline().split(" ")]
            rec["Het"] = [flt(x) for x in f.readline().split(" ")]
            yield rec
        l = f.readline()
    f.close()


def getBlocks(cfg):
    blocks = []
    gensToDo = list(cfg.futureGens)
    currBlock = list(cfg.refGens)
    blocks.append(currBlock)
    while len(gensToDo) > 0:
        if len(currBlock) == 10:
            currBlock = list(cfg.refGens)
            blocks.append(currBlock)
        currBlock.append(gensToDo[0])
        del gensToDo[0]
    return blocks
