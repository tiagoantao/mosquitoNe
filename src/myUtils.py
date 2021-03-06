try:
    import ConfigParser
except ImportError:
    import configparser as ConfigParser
import os

sysDir = "/".join(os.path.realpath(__file__).split("/")[:-2])
confDir = sysDir + "/conf/"
etcConf = sysDir + "/etc.conf"

oExpr = {
    "constant": '"samp/%f/%%d/%%d/smp-%d-%%d-%%d.txt" %% ',
    "season":   '"samp/ses/%f/%%d/%%d/smp-%d-%d-%d-%d-%%d-%%d.txt" %% ',
    "decline":  '"samp/decl/%f/%%d/%%d/smp-%d-%d-%d-%%d-%%d.txt" %% ',
}


def flt(x):
    try:
        return float(x)
    except:
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
    elif cfg.demo == "season":
        cfg.popSize = config.getint("pop", "popSize")
        cfg.seasonGen = config.getint("pop", "seasonGen")
        cfg.A = config.getint("pop", "A")
        cfg.B = config.getint("pop", "B")
        try:
            cfg.T = config.getint("pop", "T")
        except ConfigParser.NoOptionError:
            cfg.T = 12
    elif cfg.demo == "decline":
        cfg.popSize = config.getint("pop", "popSize")
        cfg.declineGen = config.getint("pop", "declineGen")
        cfg.declineSize = config.getint("pop", "declineSize")

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
    elif cfg.demo == "decline":
        rExpr = (oExpr[cfg.demo] +
                 '(numIndivs, numLoci, gen, rep)') % (
                     cfg.mutFreq, cfg.popSize, cfg.declineGen, cfg.declineSize)
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
    elif cfg.demo == "decline":
        return "samp/conc/decl-%d-%d-%f-%d-%d-%d-%d.txt" % (
            numIndivs, numLoci, cfg.mutFreq, cfg.popSize,
            cfg.declineGen, cfg.declineSize, rep)
    elif cfg.demo == "constant":
        return "samp/conc/con-%d-%d-%f-%d-%d.txt" % (
            numIndivs, numLoci, cfg.mutFreq, cfg.popSize, rep)


def getStatName(cfg, numIndivs, numLoci):
    if cfg.demo == "constant":
        return "samp/con-%d-%d-%d.txt" % (numIndivs, numLoci, cfg.popSize)
    elif cfg.demo == "season":
        return "samp/ses-%d-%d-%d-%d-%d-%d.txt" % (numIndivs, numLoci,
                                                   cfg.popSize, cfg.A,
                                                   cfg.B, cfg.T)
    elif cfg.demo == "decline":
        return "samp/decl-%d-%d-%d-%d-%d.txt" % (numIndivs, numLoci,
                                                 cfg.popSize, cfg.declineGen,
                                                 cfg.declineSize)


def getMLNE(cfg, numIndivs, numLoci):
    if cfg.demo == "constant":
        for rep in range(cfg.reps):
            for refGen in cfg.refGens:
                for futureGen in cfg.futureGens:
                    fname = "bak/mne-simple%d-%d-%d-%d-%d-%d" % (
                        cfg.popSize, rep, refGen, futureGen,
                        numIndivs, numLoci)
                    f = open(fname)
                    yield (rep, refGen, futureGen,
                           [float(x)
                            for x in f.readline().rstrip().split(' ')])
    elif cfg.demo == "decline":
        for rep in range(cfg.reps):
            for refGen in cfg.refGens:
                for futureGen in cfg.futureGens:
                    fname = "bak/mne-decl-%d-%d-%d-%d-%d-%d-%d" % (
                        cfg.popSize, cfg.declineSize, rep, refGen, futureGen,
                        numIndivs, numLoci)
                    f = open(fname)
                    yield (rep, refGen, futureGen,
                           [float(x)
                            for x in f.readline().rstrip().split(' ')])
    elif cfg.demo == "season":
        for rep in range(cfg.reps):
            for refGen in cfg.refGens:
                for futureGen in cfg.futureGens:
                    fname = "bak/mne-season-%d-%d-%d-%d-%d-%d-%d" % (
                        cfg.A, cfg.B, rep, refGen, futureGen,
                        numIndivs, numLoci)
                    f = open(fname)
                    yield (rep, refGen, futureGen,
                           [float(x)
                            for x in f.readline().rstrip().split(' ')])


def getStat(f):
    l = f.readline()  # pcrits
    l = f.readline()
    while l != "":
        rep, gen = tuple([int(x) for x in l.rstrip().split(' ')])
        rec = {}
        l = l.readline().rstrip()
        rec["rep"] = rep
        if l.find("temp") > -1:
            toks = l.split(" ")
            rec["type"] = "temp"
            rec["g1"] = int(toks[2])
            rec["g2"] = int(toks[3])
            rec["Pollak"] = [flt(x) for x in toks[4].split("#")]
            rec["Jorde/Ryman"] = [flt(x) for x in toks[5].split("#")]
            rec["Nei/Tajima"] = [flt(x) for x in toks[6].split("#")]
            yield rec
        else:
            toks = l.split(" ")
            rec["type"] = "notemp"
            rec["gen"] = gen
            rec["coanc"] = flt(toks[1])
            toks = f.readline().rstrip().split(' ')
            rec['LD'] = []
            for tok in toks[1:]:
                ld_case = tok.split('#')
                rec['LD'].append((flt(ld_case[0]), flt(ld_case[1]),
                                  flt(ld_case[2])))
            toks = f.readline().rstrip().split(' ')
            rec["Het"] = [flt(x) for x in toks[1]]
            yield rec
        l = f.readline()
    f.close()
