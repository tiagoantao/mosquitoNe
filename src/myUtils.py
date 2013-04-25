import ConfigParser
import os

sysDir = "/".join(os.path.realpath(__file__).split("/")[:-2])
confDir = sysDir + "/conf/"
etcConf = sysDir + "/etc.conf"

oExpr = {
    "constant": '"samp/%f/%%d/%%d/smp-%d-%%d-%%d.txt" %% ',
    "season":   '"samp/ses/%f/%%d/%%d/smp-%d-%d-%d-%%d-%%d.txt" %% ',
}


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
                     cfg.mutFreq, cfg.popSize, cfg.A, cfg.B)
    else:
        rExpr = (oExpr[cfg.demo] +
                 '(numIndivs, numLoci, gen, rep)') % (cfg.mutFreq, cfg.popSize)
    if doeval:
        return eval(rExpr)
    else:
        return rExpr


def getConc(cfg, numIndivs, numLoci, rep):
    if cfg.demo == "season":
        return "samp/ses-%d-%d-%f-%d-%d-%d-%d.txt" % (numIndivs, numLoci,
                                                      cfg.mutFreq, cfg.popSize,
                                                      cfg.A, cfg.B, rep)
    else:
        return "samp/con-%d-%d-%f-%d-%d.txt" % (numIndivs, numLoci,
                                                cfg.mutFreq, cfg.popSize, rep)


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
