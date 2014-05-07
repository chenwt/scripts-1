# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <rawcell>

# This is to use greedy methods to identify a subset of ceRNA regulator * mutated samples for one target cancer gene. 
# This subset should optimized the correlation between sum-up ceRNA driver expression and target expression.
# steps involved:
# 
# 
# 'reduceColumn method'
# -1 prepareData
# -2 calculated fullMatrix correlation corr_0, fisher transformed 
# -3 start, set corr_0 = 0, set one sample inactive
# -4 calculate corr_k, fisher transformed, 
# -4.5 if smpss_k != smpss_k-1: corr_k = fisher.transfor(corr_k)
# -5 if corr_k >= corr_k-1: deactivate mutation , otherwise keep it 
# -6 repeat step4-5,
#    
# 
# 'Cell method'
# -1 prepareData
# -2 calculated fullMatrix correlation corr_0, fisher transformed 
# -3 set corr_k-1 = corr_prev, flip one mutation from 1 to 0(deactive it)
# -4 calculate corr_k, fisher transformed, 
# -5 if corr_k - corr_k-1 > tolenrance: deactivate mutation , otherwise filp it back to 1 (active mutation) 
# -6 repeat step4-5 until all mutations have been visited.
# 
# 'select optimal tolenrance'
# The idea: the more stringent tol produce less mutation.
# select the tolerance which produe maxi correlation, with significant p-value compare to full matrix and permut matrix
# 
# 'select optimal over all random initiation' 
# The idea: random initiation will proceduce slight different result, so try to approximate to global optmization by randomization.
# -1 Ran multiply random initization ( n = 100), select the one with optimized correlation, with qualify p-values
# -2 count the number of mutations being selected, output the final mutation


keyRegSumfile="/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/runMay5/keyRegSumfile.cancergene.05052014"
expfile="/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_tumor_Mar-21-2014.matrix_Mar-26-2014.voomNormed.matrix"
gslistfile="//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/gslist/gslist_CnvMethSomFree.10smapMore.deg_20140430.txt.hasReg.list"
mutfile="/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/tcgal2som/brca_somlevel2_byGene.matrix.inExpSmp.matrix.nonzero"
output="/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/test/test_optCorr_05062014"

import  sys,getopt
from    collections     import defaultdict
from    parseKeyRegFile import parseKeyRegFile
from    collections     import Counter, Sequence
from    parseGslistFile import parseGslistFile
import  numpy           as np
import  subprocess
import  os.path


def formatSampleName(code19):
    if len(code19) >11:
        return code19[5:16].replace("-",".")
    else :
        return code19.replace("-", ".")


def loadTarReg(keyRegSumfile):
    resD = defaultdict(list)
    with open(keyRegSumfile) as f:
        line = f.readline()
        while line:
            crt_t, ctr_r = line.strip().split("\t")
            resD[crt_t] = ctr_r.split(";")
            line = f.readline()
    return(resD)


def getMutExpSmp(expfile, mutfile):
    resL = []
    with open(expfile) as f:
        resL = f.readline().strip().split("\t")
        resL = map(formatSampleName, resL)
    with open(mutfile) as f:
        allMutSmp = f.readline().strip().split("\t")
        resL = [a for a in map(formatSampleName, allMutSmp) if a in resL]
    return resL


def loadTarIntSmp(gslistfile):
    resD = defaultdict(list)
    with open(gslistfile) as f:
        line = f.readline()
        while line:
            crt_t, crt_s = line.strip().split("\t")
            resD[crt_t] = crt_s.split(";")
            line = f.readline()
    return(resD)


def loadExp(expfile, smpsL):
    resD = defaultdict(list)
    with open(expfile) as f:
        allExpSmp = map(formatSampleName, f.readline().strip().split("\t"))
        smpIndx = [id for (id, v) in enumerate(allExpSmp) if v in smpsL]
        resD['gene'] = map(allExpSmp.__getitem__, smpIndx)

        line = f.readline()
        while line:
            crt_g, crt_e = line.strip().split("\t",1)
            temp = map(float,crt_e.split("\t"))
            resD[crt_g] = map(temp.__getitem__, smpIndx)
            line = f.readline()
    return resD

def loadMut(mutfile, smpsL):
    resD = defaultdict(list)
    with open(mutfile) as f:
        _, allMutSmp = f.readline().strip().split("\t",1)
        allMutSmp = map(formatSampleName, allMutSmp.split("\t"))
        smpIndx = [id for (id, v) in enumerate(allMutSmp) if v in smpsL]
        resD['gene'] = map(allMutSmp.__getitem__, smpIndx)
        line = f.readline()
        while line:
            crt_g, crt_v = line.strip().split("\t",1)
            temp =  map(int,crt_v.split("\t"))
            resD[crt_g] = map(temp.__getitem__, smpIndx)
            line = f.readline()
    return resD

tarRegD = loadTarReg(keyRegSumfile)
tarIntSmpD = loadTarIntSmp(gslistfile)
mutExpsmpL = getMutExpSmp(expfile, mutfile)
expD = loadExp(expfile, mutExpsmpL)
mutD = loadMut(mutfile, mutExpsmpL)

for tgene in tarRegD.keys(): 
    # def startOpt4Gene(tgene, tarIntSmpD, tarRegD, mutD, expD, output):
    tIntSmp = tarIntSmpD[tgene]
    allRegsL = tarRegD[tgene]
    intMutSmpIdL = [id for (id, s) in enumerate(mutD['gene']) if s in tIntSmp]
    intExpSmpIdL = [id for (id, s) in enumerate(expD['gene']) if s in tIntSmp]
    
    regMutD = {k:map(v.__getitem__, intMutSmpIdL) for (k,v) in mutD.items() if k in allRegsL}
    regExpD = {k:map(v.__getitem__, intExpSmpIdL) for (k,v) in expD.items() if k in allRegsL}
    
    expMutRegL = set(regExpD.keys()).intersection(set(regMutD.keys()))
    outTempMut = output + "_" + tgene + "_regMut.temp"
    outTempExp = output + "_" + tgene + "_exp.temp"
    print outTempMut 
    print outTempExp 
    outTempMutH = open(outTempMut,'w')
    outTempExpH = open(outTempExp,'w')
    
    outTempMutH.write('gene\t'+"\t".join( map(mutD['gene'].__getitem__, intMutSmpIdL)) + "\n")
    for k,v in regMutD.items():
        if k in expMutRegL: 
            outTempMutH.write(k + "\t" + "\t".join(map(str,v)) + "\n")
    
    outTempExpH.write('gene\t'+"\t".join(map(expD['gene'].__getitem__, intExpSmpIdL)) + "\n")
    
    ## target expession in the first row
    outTempExpH.write(tgene + "\t" + "\t".join(map(str,map(expD[tgene].__getitem__, intExpSmpIdL))) +"\n")
    for k,v in regExpD.items():
        if k in expMutRegL:
            outTempExpH.write(k + "\t" + "\t".join(map(str,v)) + "\n")  
    outTempMutH.close()
    outTempExpH.close()
    
    outTemp = output + "_" + tgene + ".tsv"
    #     cmd = "/ifs/home/c2b2/ac_lab/jh3283/tools/R/R-3-02/bin/Rscript \
    #     /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step2-2_regKeyRegulators_v5.r\
    #     --vanilla --exp " +  outTempExp + " --mut " + outTempMut + " --output " + outTemp 
    cmd = "/ifs/home/c2b2/ac_lab/jh3283/tools/R/R-3-02/bin/Rscript /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step3-4_greedyOptCorr.r \
    --vanilla --exp " +  outTempExp + " --mut " + outTempMut + " --output " + outTemp 
    if not os.path.isfile(outTemp) : 
        subprocess.Popen(cmd, shell = True)
    else:
        print outTemp + " already exist, skip it " 
#     rtncode  = subprocess.call(cmd, shell = True)
#     if rtncode != 0 :
#         print "Error in regression!"
#         sys.exit()
#     print tgene + " Done"


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

def format2Dlist(l2d):
    out = []
    for i in l2d:
        out.append("["+",".join(map(str,i)) + "]")
    return ";".join(out)


