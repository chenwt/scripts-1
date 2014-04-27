import numpy as np
import  sys,getopt
from    collections     import defaultdict
from    parseKeyRegFile import parseKeyRegFile
from    collections     import Counter, Sequence
from    parseGslistFile import parseGslistFile


class MutSet():
    def __init__(self, gene, mutSample, zscore):
        self.gene = gene
        self.mutSample = mutSample
        self.mutZscore = dict(zip(mutSample,zscore))
    def __str__(self):
        return "(%s, %s)" % (self.gene, self.mutSample)
    def __expr__(self):
        return "(%s, %i mutSmp)" % (self.gene, len(self.mutSample))
    def update(self, unionSample):
        for smp in self.mutSample:
            if not smp in unionSample:
                self.mutSample = self.mutSample.remove(smp)
                self.mutZscore = self.mutZscore.remove(smp)
            else:
                pass

def findWegitedMin(S, R, wtype = 'mean'):
    ''' 
    wtype is to select the method to use weight,
    total: the summation of all mutation zscores;
    mean: the mean of all mutation zscores;
    max: the maximization of mutation zscores for each gene
    '''
    ## get the minimum cost set
    minCost = 99999.0
    minElement = -1
    minSet = ''
    minSname = ''
    for i, s in enumerate(S):
        sname_i = s.gene; ss = s.mutSample; sw = s.mutZscore
        ss_i  = set(R).intersection(set(ss))
        if len(ss_i) == 0:
            continue
        sw_i = [sw[a] for a in ss_i ]
        if wtype == 'total':
            cost = 1/reduce(lambda x, y: x+y , sw_i)
        elif wtype == 'mean':
            cost = len(sw_i)/sum(sw_i)
        elif wtype == 'max':
            cost = 1/ reduce(lambda x, y: max(x,y), sw_i)
        else:
            print "ERRor wtype: use default mean; other option total, max or mean"
            cost = len(sw_i)/sum(sw_i)
        if cost < minCost:
            minCost = cost
            minElement = i
            minSname  = sname_i
            minSet  = ss_i 
    return minSname, minSet, minCost
def wgsc(S, U, alpha = 0.8, tol = 0.001, wtype = 'mean'):
    R = U
    C = []
    G = []
    costs = []
    while len(R) != 0:
       g_i, S_i, cost = findWegitedMin(S, R, wtype = wtype)
       C.append(list(S_i))
       G.append(g_i)
       R = list(set(R).difference(set(S_i)))
       costs.append(cost)
       if len(R) < int((1 - alpha) * len(U)):
            break
    return G, C, costs

def __test__():
    from collections import defaultdict
    seqSet  = {'G1':['S2','S4','S6'],
                'G2':['S1','S3'],
                'G3':['S1'],
                'G4':['S1'],
                'G5':['S5'],
                'G6':['S3']}
    seq     = ['S1', 'S2', 'S3', 'S4', 'S5','S6']
    weightSet = {'G1':[1.0,0.5,1.5],
                 'G2':[2.0, 2.5],
                 'G3':[2.3],
                 'G4':[1.2],
                 'G5':[2.5],
                 'G6':[3.0]}
    setObjL = []
    for sk, ss in seqSet.items():
        setObjL.append(MutSet(sk,ss,weightSet[sk]))
    geneL, smpL, costL = wgsc(setObjL, seq, wtype = "mean")
    geneL, smpL, costL = wgsc(setObjL, seq, wtype = "total")
    geneL, smpL, costL =  wgsc(setObjL, seq, wtype = "max")

####-------compute the zscore matrix for each mutation
import numpy as np
# def myZscore(a,b):
#     bArr = np.array(b)
#     m  = np.mean(bArr)
#     sd = np.std(bArr)
#     return (np.array(a) - m)/sd
def myZscore(a,b):
    return abs((np.array(a) - b[0])/b[1])

def formatSampleName(code19):
    if len(code19) >11:
        return code19[5:16].replace("-",".")
    else :
        return code19.replace("-", ".")


def loadNormExpfile(filename):
    expD = defaultdict(list)
    with open(filename) as f:
        samples = f.readline().strip().split("\t")
        line = f.readline()
        while line:
            gCrt, valCrt = line.strip().split("\t",1)
            valCrt = np.array(map(float, valCrt.split("\t")))
            expD[gCrt] = [np.mean(valCrt), np.std(valCrt)]
            line = f.readline()
    return expD
  
def loadExpfile(filename, expND):
    expD = defaultdict(list)
    with open(filename) as f:
        expD['samples'] = f.readline().strip().split("\t")
        line = f.readline()
        while line:
            gCrt, valCrt = line.strip().split("\t",1)
            try: 
                expD[gCrt] = map(lambda x:\
                                 myZscore(float(x),expND[gCrt]),\
                                 valCrt.split("\t"))
            except:
                pass 
            line = f.readline()
    return expD

def loadMutInfo(mutfile, zscoreD):
    ''' 
    load all mutated ceRNA driver, 
    and return all mutated ceRNA driver's mutated sample, and zscores
    '''
    mutD = defaultdict(list)
    mutZscoreD = defaultdict(list)
    cnt  = 0 
    with open(mutfile) as f:
        gene, samples = f.readline().strip().split("\t",1)
        samples =  map(formatSampleName, samples.split("\t"))
        mutD['samples'] = samples
        line = f.readline()
        while line:
            cnt = cnt + 1
            gene, vals = line.strip().split("\t",1)
            mutIdx = [id for (id,m) in enumerate(vals.split("\t")) if m != "0"]
            mutSmp = map(samples.__getitem__, mutIdx)
            mutSmpInExpSmpID = []; mutSmpInExpSmp =  []
            for (id, a) in enumerate(zscoreD['samples']) :
                if a in mutSmp:
                    mutSmpInExpSmpID.append(id)
                    mutSmpInExpSmp.append(a)
            mutZscoreD[gene] = map(zscoreD[gene].__getitem__, mutSmpInExpSmpID)
            mutD[gene] = mutSmpInExpSmp
            line = f.readline()
    print " input target genes:\t",cnt 
    return mutD, mutZscoreD

def prepareDataWGSC(mutfile, gslistfile, keygenefile, pvalCut = 0.01 ):
    '''
    given, mutation dict, zscore dict, 
    intact sample list for each cancer genen, 
    file from group lasso result; 
    prepare mutation data, zscore data for each cancer gene
    using MutSet objects
    '''
    
    tgeneSum, regsSum   = parseKeyRegFile(keygenefile, pvalCut)
    if not regsSum or not tgeneSum:
        return tgeneSum, ''
    reglist         = regsSum
    regMutObjL      = [] 
    tgene           = tgeneSum[0]
    gintsmplist     = parseGslistFile(tgene, gslistfile)
   
    ##---check whether mutated samples are intact
    cnt = 0 
    regMutAllSmp = []
    mutRegs = []
    for gene in reglist:
        cnt = cnt + 1
        crtMut = mutD[gene]; crtMutZscore = mutZscoreD[gene]
        if crtMut:
            mutRegs.append(gene)
        for idx, smp in enumerate(mutD[gene]) :
            if smp not in gintsmplist:
                crtMut.remove(smp)
                del crtMutZscore
            else:
                regMutAllSmp.append(smp)
                pass
        if crtMut:
            regMutObjL.append(MutSet(gene, crtMut, crtMutZscore))
    tempCnter = Counter(regMutAllSmp)
    regMutAllSmp, regMutAllSmpLoad = tempCnter.keys(), tempCnter.values()
    outInfo = {}
    outInfo['gintSmp']          = gintsmplist       ## all gint sample 
    outInfo['mutGintSmp']       = regMutAllSmp     ## mutated gint sample 
    outInfo['mutGintSmpLoad']   = regMutAllSmpLoad  ## mutated gint sam  ple's mutation frequency for each mutation
    outInfo['tgene']            = tgene     ## target gene
    outInfo['allRegs']          = reglist  ## all regulators
    outInfo['mutRegs']          = mutRegs## all mutated regulators
    outInfo['intMutRegs']       = map(lambda x:x.gene, regMutObjL)
    return  {tgene:[regMutObjL, regMutAllSmp]}, outInfo


def __test__():
    expTumorM="/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_tumor_Mar-21-2014.matrix_Mar-26-2014.voomNormed.matrix"
    expNormalM="/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_normal_Mar-21-2014.matrix_Mar-26-2014.voomNormed.matrix"
    expnD  = loadNormExpfile(expNormalM)
    zscoreD  = loadExpfile(expTumorM, expnD)
    
    mutfile = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/sigMut/step2_mutKeyReg/kegRegs_Apr-18-2014.summary.driverRegs.list.uniq.mut.matrix"
    
    mutD, mutZscoreD = loadMutInfo(mutfile, zscoreD) 
    
    gslistfile="/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/gslist/gslist_Mar-24-2014_CnvMethSomFree.10smapMore.deg_20140325.txt.10more.hasReg.list"
    keygenefile="/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/candiReg/run-Apr-1-2014/data/BCL9_candidateRegs_Mar-31-2014.txt"
    pvalCut=0.01
    
    tRegMutObjDict, info = prepareDataWGSC(mutfile, gslistfile, keygenefile)
    targ = tRegMutObjDict.keys()[0]
    print wgsc(tRegMutObjDict[targ][0],tRegMutObjDict[targ][1] , wtype = "mean")
    print wgsc(tRegMutObjDict[targ][0],tRegMutObjDict[targ][1] , wtype = "total")
    print wgsc(tRegMutObjDict[targ][0],tRegMutObjDict[targ][1] , wtype = "max")


