# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <rawcell>

# To prepare data for running regression:
# 1. TF target CNV
# 2. TF expression  
# 3. TF target mutation profile from TF greedy

# <codecell>

import  sys,getopt
from    collections     import defaultdict, Counter, Sequence
import  subprocess
import  os.path
import  pandas as pd
from    pandas import DataFrame, Series
import  re

# <codecell>

tfGdResfile = '/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/reg/tfTar/sigTgTFSmp.txt'
smpOrderfile = '/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/reg/tfTar/smp.order'

cnvfile = '/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/cnv/brca_cnv_l3_tumor_Mar-23-2014.matrix.uniq.matrix'
expfile = '/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_tumor_Mar-21-2014.matrix_Mar-26-2014.voomNormed.matrix'

output = '/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/reg/tfTar/test/sigTgTFSmp.TFmutProfile'
outputDir = '/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/reg/tfTar/runAug9/'

# <codecell>

def name2binSL(tfsmps):
    tf, smps = tfsmps[1:-1].split(":")
    smps = smps.split(",")
    binSmpsl = Series([0] * len(smpOrder), index=smpOrder)
    binSmpsl.loc[smps] = 1
    return tf, binSmpsl

# <codecell>

def rec2TfSmpDF( tfSmpAll):
    val = tfSmpAll
    ValLst = val.split(";")
    tempDt = dict(map(name2binSL, ValLst))
    tfs = tempDt.keys()
    smpSLs = tempDt.values()
    del tempDt
    smpDF = pd.concat(smpSLs, axis=1, keys = tfs)
    return smpDF

# <codecell>

def formatSampleName(code19):
    if len(code19) >11:
        return code19[5:16].replace("-",".")
    else :
        return code19.replace("-", ".")

# <codecell>

def loadExp(expfile, smpsL):
    resD = defaultdict(list)
    with open(expfile) as f:
        allExpSmp = map(formatSampleName, f.readline().strip().split("\t"))
        line = f.readline()
        while line:
            crt_g, crt_e = line.strip().split("\t",1)
            temp = map(float,crt_e.split("\t"))
            resD[crt_g] = temp
            line = f.readline()
    resDF = DataFrame(resD, index=allExpSmp).T[smpsL]
    return resDF

# <codecell>

def loadCNV(cnvfile, smpsL):
    resD = defaultdict(list)
    with open(cnvfile) as f:
        allSmp = map(formatSampleName, f.readline().strip().split("\t")[1:])
        commSmp = list(set(smpsL).intersection(set(allSmp)))
#         print len(allSmp), len(commSmp), len(smpsL)
        
        line = f.readline()
        while line:
            crt_g, crt_e = line.strip().split("\t",1)
            temp = map(float,crt_e.split("\t"))
            resD[crt_g] = temp
            line = f.readline()
    resDF = DataFrame(resD, index=allSmp).T[commSmp]
    return resDF

# <codecell>

def prepFile(tfGdResfile, commSmp, expDF, cnvDF):
    with (open(tfGdResfile)) as f:
        gCNV = list(cnvDF.index)
        
        line = f.readline()
        while line:
            tg_crt, tfSmpAll = line.strip().split()
            
            if not tg_crt in gCNV :
                line = f.readline()
                continue
    
            tfSmpDF = rec2TfSmpDF(tfSmpAll).T            
            commSmp = list(set(tfSmpDF.columns).intersection(set(commSmp)))
                
            tmp = pd.concat([expDF.loc[tg_crt,commSmp], cnvDF.loc[tg_crt,commSmp]], axis=1).T
            tmp.index = ['exp','cnv']
            
            outDF = pd.concat([tmp, tfSmpDF[commSmp]])
            
    #             tfSmpDF.columns = map(lambda x: tg_crt + "_" + x, list(tfSmpDF.columns))
            outDF.to_csv(outputDir + "/" + tg_crt + ".tgTFreg.input" ,\
                             sep="\t",quoting=False, mode = 'w')
            line = f.readline()

# <codecell>

# tg = 'IFFO2'

# tmp.iloc[:,:5]

# <codecell>

### main

#get ordered sample name
smpOrder = []
with (open(smpOrderfile)) as f:
    smpOrder = f.readline().strip().split("\t")
##load data
expDF = loadExp(expfile, smpOrder)
cnvDF = loadCNV(cnvfile, smpOrder)
commSmp = list(set(expDF.columns).intersection(set(cnvDF.columns)))

##output
prepFile(tfGdResfile, commSmp, expDF, cnvDF)

# <codecell>


# <codecell>


# <codecell>


