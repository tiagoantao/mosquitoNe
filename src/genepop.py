import random
import exceptions

def save(pop, numIndivs, numLoci, oExpr, ploidy):
    fName = oExpr
    try:
        f = open(fName, "w")
    except exceptions.IOError:
        raise exceptions.IOError, "Can not open file " + fName + " to write."
    #
    # file is opened.
    np = pop.numSubPop()

    loci = range(numLoci)
    nl = len(loci)

    nd = 3
    # write the first line
    f.write( 'Sim %d %d %d %d\n' % (np, numIndivs, numLoci, ploidy) )
    # following lines with loci name.
    for li in range(len(loci)):
        lname = pop.locusName(loci[li])
        if lname == "":
            lname = "locus_" + str(li+1)
        f.write( lname +"\n");
    totNumLoci = pop.totNumLoci()
    for sp in range(0, pop.numSubPop()):
        # genotype of subpopulation sp, individuals are
        # rearranged in perfect order
        # here is where the age classes infoField should be
        f.write('Pop\n')
        gt = pop.genotype(sp)
        indivs = random.sample(range(pop.subPopSize(sp)), numIndivs)
        for ind in indivs:
            i = pop.individual(ind, sp)
            f.write("%d(%d), " % (sp,ind))
            if ploidy==2:
                p1 = 2*totNumLoci*ind          # begining of first hemo copy
                p2 = 2*totNumLoci*ind + totNumLoci     # second
            else:
                p1 = totNumLoci*ind          # begining of first hemo copy
            for al in loci: # allele
                ale1 = gt[p1+al]
                if ploidy==2:
                    ale2 = gt[p2+al]
                    f.write('%%0%dd%%0%dd ' % (nd, nd) % (ale1+1, ale2+1))
                else:
                    f.write('%%0%dd ' % (nd) % (ale1+1))
            f.write( "\n")
    f.close()
