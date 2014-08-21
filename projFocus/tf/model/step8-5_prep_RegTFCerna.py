# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <rawcell>

# To prepare data for running regression:

# <codecell>

sys.modules[__name__].__dict__.clear()
import  sys,getopt
from    collections     import defaultdict, Counter, Sequence
import  subprocess
import  os.path
import  pandas as pd
from    pandas import DataFrame, Series
import  re
import  sys
from    math import copysign

# <codecell>

def sign(x): 
    if x > 0.0 : 
        return 1
    elif x == 0.0 or x == 0: 
        return 0 
    else: 
        return -1 

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

def maxAbs(x,y):
    sx = sign(x); sy = sign(y)
    x = abs(x); y = abs(y)
    if x >=y:
        return(sx * x) 
    else:
        return(sy * y)
    
def loadCNV(cnvfile):
    resD = defaultdict(list)
    with open(cnvfile) as f:
        allSmp = map(formatSampleName, f.readline().strip().split("\t")[1:])
#         print len(allSmp), len(commSmp), len(smpsL)
        line = f.readline()
        while line:
            crt_g, crt_e = line.strip().split("\t",1)
            temp = map(float,crt_e.split("\t"))
            if crt_g in resD.keys(): 
                resD[crt_g] = map(lambda x: maxAbs(x[0],x[1]), zip(resD[crt_g],temp))
            else:
                resD[crt_g] = temp
                
            line = f.readline()
    resDF = DataFrame(resD, index=allSmp).T
    return resDF

# <codecell>

# def writeInputForRegCerna(cTar_crt):
    
if not cTar_crt in geneExp :
#         print "no cTar expression"
    return 0

cRegs_crt = cTarsRegDt[cTar_crt]
exp_cTar_crt = expDF.loc[[cTar_crt] + list(geneExp.intersection(set(cRegs_crt))), smpLst].T;
exp_cTar_crt.columns = ['cTarExp'] + map(lambda x: "cRegExp_" + x, cRegs_crt)

#     if cTar_crt in geneCnv:
#         cnv_cTar_crt = cnvDF.loc[[cTar_crt],smpLst].T; cnv_cTar_crt.columns = ['cTarCNV']
#     else:
#         print "no selfCNV"

# cRegs_crt = list(set(geneCnv).intersection(set(cRegs_crt)))
cnv_cReg_crt = cnvDF.loc[set(geneCnv).intersection(set(cRegs_crt)), smpLst].T; 
cnv_cReg_crt.columns = map(lambda x: 'cRegCNV_' + x, cnv_cReg_crt.columns)

### collapsed TF profiles
# # # cRegs_crt = list(set(geneTFmut).intersection(set(cRegs_crt)))
mut_cRegTF_crt = mutOptMergeDF.loc[cRegs_crt,smpLst].T;
mut_cRegTF_crt.columns = map(lambda x: 'cRegTFmut_' + x, mut_cRegTF_crt.columns)



### prepare cReg TF's mutation matrix

cRegs_crt = cTarsRegDt[cTar_crt]; 
out_tfMutOptDF = pd.concat([exp_cTar_crt, cnv_cReg_crt], axis = 1)
cntCol1 = len(list(out_tfMutOptDF.columns))
for tmp_rg in list(set(cRegs_crt).intersection(set(cRegsTfDt.keys()))): 
    cReg_TFs_crt = cRegsTfDt[tmp_rg]
    
    # print cTar_crt, cRegs_crt, tmp_rg, cReg_TFs_crt
    
    
    tmp_tfMutAct = mutActDF.loc[cReg_TFs_crt,smpLst]
    tmp_cRegOptMut = mutOptMergeDF.loc[tmp_rg, smpLst]
    
#         mutOptMergeDF.loc[[tmp_rg], smpLst].to_csv(output + ".optIntMut.temp", sep="\t")
    
#         tmp_tfMutAct.to_csv(output + ".actMut.temp", sep="\t")
    noMutSmp = set(tmp_tfMutAct.columns).difference(set(tmp_cRegOptMut[tmp_cRegOptMut > 0].index))
    tmp_tfMutAct.loc[:,noMutSmp] = 0
    tmp_tfMutAct.index =  map(lambda x: 'cRegTFmut_' + x + "_" + tmp_rg, tmp_tfMutAct.index)
    out_tfMutOptDF = pd.concat([out_tfMutOptDF, tmp_tfMutAct.T], axis=1)

cntCol2 = len(list(out_tfMutOptDF.columns))

#     if cntCol2 > cntCol1 : 
#         out_tfMutOptDF.to_csv(output + cTar_crt, sep="\t")
#         return 1
#     else :
# #         print cTar_crt + " no TFmut"
#         return 0 
    

# <codecell>


# <codecell>

gfile = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/reg/tfCerna/reg_tfTarMut_runAug13.txt.significant"

cnvfile = '/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/cnv/brca_cnv_l3_tumor_Mar-23-2014.matrix.uniq.matrix'
expfile = '/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_tumor_Mar-21-2014.matrix_Mar-26-2014.voomNormed.matrix'

mutActfile="/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/tcgal2som/brca_somlevel2_byGene.matrix"
mutOptfileCollapse="/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/tfMut/runJuly27/summary/all_targets_mutSampleVector"

cernafile = '/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/data/networks/brca_ceRNA_network.txt'
cTarcRegGrplasfile = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/candiReg/runGATA3/summary/tf_aracne_gata3_0.01_0"


output = '/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/reg/tfCerna/data/input_'
outputDir = '/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/reg/tfCerna/data/'

# <codecell>


# <codecell>

### main-------
### given cerna regulators, extract all it's target, and create tar-Reg dictionary
expDF = pd.read_csv(expfile, delimiter="\t", index_col=0); expDF.columns = map(formatSampleName, expDF.columns)
cnvDF = loadCNV(cnvfile)

cRegsTfDt = defaultdict(list)
with(open(gfile)) as f:
    line = f.readline()
    while line:
        cReg_crt, tfs_crt = line.strip().split("\t")
        tfs_crt = tfs_crt.split(",")
        if len(tfs_crt) > 0 :
            cRegsTfDt[cReg_crt].extend(tfs_crt)
        line = f.readline()

mutActDF = pd.read_csv(mutActfile, index_col=0, delimiter="\t")
mutActDF.columns = map(formatSampleName, list(mutActDF.columns))

mutOptMergeDF = pd.read_csv(mutOptfileCollapse,delimiter="\t",index_col = 0)
mutOptMergeDF.columns = map(formatSampleName, mutOptMergeDF.columns)

cRegsLst = list(mutOptMergeDF.index)

geneCnv = list(cnvDF.index); geneTFmut = list(mutOptMergeDF.index)


### cernaTar---ceranReg----
### From extract genes from 
cTarsRegDt = defaultdict(list)
with(open(cernafile)) as f:
    line = f.readline()
    while line:
        r1, r2, _ = line.strip().split("\t", 2)
        if not ( r1 in cRegsLst or r2 in cRegsLst):
            line = f.readline()
            continue
        elif r1 in cRegsLst:
            cTarsRegDt[r2].append(r1)
        elif r2 in cRegsLst:
            cTarsRegDt[r1].append(r2)
        line = f.readline()
  
cTarsLst = list(cTarsRegDt.keys()) 
### ---------


smpLst = list(set(list(expDF.columns)).intersection(set(list(mutOptMergeDF.columns))).intersection(set(list(cnvDF.columns))))

# <codecell>

geneExp = set(expDF.index)
geneActMut = set(mutActDF.index)

def writeInputForRegCerna(cTar_crt):
    
    if not cTar_crt in geneExp :
        return 0
    cRegs_crt = cTarsRegDt[cTar_crt]
    exp_cTar_crt = expDF.loc[[cTar_crt] + list(geneExp.intersection(set(cRegs_crt))), smpLst].T;
    exp_cTar_crt.columns = ['cTarExp'] + map(lambda x: "cRegExp_" + x, cRegs_crt)
    
#     if cTar_crt in geneCnv:
#         cnv_cTar_crt = cnvDF.loc[[cTar_crt],smpLst].T; cnv_cTar_crt.columns = ['cTarCNV']
#     else:
#         print "no selfCNV"

    cRegs_crt = list(set(geneCnv).intersection(set(cRegs_crt)))
    cnv_cReg_crt = cnvDF.loc[set(geneCnv).intersection(set(cRegs_crt)), smpLst].T; 
    cnv_cReg_crt.columns = map(lambda x: 'cRegCNV_' + x, cnv_cReg_crt.columns)
    
    ### collapsed TF profiles
    # # # cRegs_crt = list(set(geneTFmut).intersection(set(cRegs_crt)))
    mut_cRegTF_crt = mutOptMergeDF.loc[cRegs_crt,smpLst].T;
    mut_cRegTF_crt.columns = map(lambda x: 'cRegTFmut_' + x, mut_cRegTF_crt.columns)
    
    
    
    ### prepare cReg TF's mutation matrix
    cRegs_crt = cTarsRegDt[cTar_crt]; 
    out_tfMutOptDF = pd.concat([exp_cTar_crt, cnv_cReg_crt], axis = 1)
    cntCol1 = len(list(out_tfMutOptDF.columns))
    
    for tmp_rg in list(set(cRegs_crt).intersection(set(cRegsTfDt.keys()))): 
        ## for each cRegulator, extract all it's TF optmized mutation profile
        cReg_TFs_crt = cRegsTfDt[tmp_rg]
        
        tmp_tfMutAct = mutActDF.loc[geneActMut.intersection(set(cReg_TFs_crt)),smpLst]
        tmp_cRegOptMut = mutOptMergeDF.loc[tmp_rg, smpLst]
        
        
        noMutSmp = set(tmp_tfMutAct.columns).difference(set(tmp_cRegOptMut[tmp_cRegOptMut > 0].index))
        tmp_tfMutAct.loc[:,noMutSmp] = 0
        tmp_tfMutAct.index =  map(lambda x: 'cRegTFmut_' + tmp_rg  + "_" +  x, tmp_tfMutAct.index)
        out_tfMutOptDF = pd.concat([out_tfMutOptDF, tmp_tfMutAct.T], axis=1)
    
    cntCol2 = len(list(out_tfMutOptDF.columns))
    
    if cntCol2 > cntCol1 : 
        out_tfMutOptDF.to_csv(output + cTar_crt, sep="\t")
        return 1
    else :
        return 0 
    
    

# <codecell>

# len(list(cnvDF.index))
# map(writeInputForRegCerna, cTarsLst)

cntFail = 0; cntSuc = 0 
for cTar_crt in cTarsLst:
    reOut = writeInputForRegCerna(cTar_crt)
    if reOut == 0:
        cntFail = cntFail + 1
    else:
        cntSuc = cntSuc + 1 
print "fail: " + str(cntFail) + " output number " + str(cntSuc)

# <codecell>


# <codecell>


# <codecell>


# <codecell>


# <codecell>


# <codecell>


