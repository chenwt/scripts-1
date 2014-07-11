# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <rawcell>

# This is to use greedy methods to identify a subset of ceRNA regulator * mutated samples for one target cancer gene. 
# This subset should optimized the correlation between sum-up ceRNA driver expression and target expression.
# steps involved:
# 
# Python code:
# -1 extract exp, mutMatrix, tarexp data
# -2 run greedy(by submitting a job to cluster, with 1000 random starting), Call R code.
# -3 permuat mutMatrix, and do step 2
# -4 do step3 1000 times
# -5 fire a job to moniter the finishing status of the submitted job
# -5 integrate data for step1, step4
# -6 generate summarized result. 
# -7 clean folder
# 
# 'Cell method'
# -1 prepareData
# -2 calculated fullMatrix correlation corr_0, fisher transformed 
# -3 set corr_k-1 = corr_prev, flip one mutation from 1 to 0(deactive it)
# -4 calculate corr_k, fisher transformed, 
# -5 if corr_k - corr_k-1 > tolenrance: deactivate mutation , otherwise filp it back to 1 (active mutation) 
# -6 repeat step4-5 until all mutations have been visited.
# 
# 'select optimal tolenrance'- flex, combine both p-value
# The idea: the more stringent tol produce less mutation.
# select the tolerance which produe maxi correlation, with significant p-value compare to full matrix and permut matrix
# 
# 'select optimal over all random initiation' 
# The idea: random initiation will proceduce slight different result, so try to approximate to global optmization by randomization.
# -1 Ran multiply random initization ( n = 100), select the one with optimized correlation, with qualify p-values
# -2 count the number of mutations being selected, output the final mutation

# <codecell>

cnvfile = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/cnv/brca_cnv_l3_tumor_Mar-23-2014.matrix.uniq.matrix"

# <codecell>

keyRegSumfile="/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/runMay5/keyRegSummary_donejob_05172014_0.01"
expfile="/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_tumor_Mar-21-2014.matrix_Mar-26-2014.voomNormed.matrix"
gslistfile="/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/gslist/gslist_CnvMethSomFree.10smapMore.deg_20140430.txt.hasReg.list"
mutfile="/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/tcgal2som/brca_somlevel2_byGene.matrix.inExpSmp.matrix.nonzero"
output="/Users/jh3283/HOME/DATA/projFocus/result/05012014/sigMut/test/test_optCorr_wCNV"
figd="/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/report/Jul2014/fig/"

# <codecell>

import  sys,getopt
from    collections     import defaultdict
from    parseKeyRegFile import parseKeyRegFile
from    collections     import Counter, Sequence
from    parseGslistFile import parseGslistFile
import  numpy           as np
import  subprocess

# <codecell>

def formatSampleName(code19):
    if len(code19) >16:
        return code19[5:16].replace("-",".")
    else :
        return code19.replace("-", ".")

# <codecell>

def loadTarReg(keyRegSumfile):
    resD = defaultdict(list)
    with open(keyRegSumfile) as f:
        line = f.readline()
        while line:
            crt_t, ctr_r = line.strip().split("\t")
            resD[crt_t] = ctr_r.split(";")
            line = f.readline()
    return(resD)

# <codecell>

def getMutExpSmp(expfile, mutfile):
    resL = []
    with open(expfile) as f:
        resL = f.readline().strip().split("\t")
        resL = map(formatSampleName, resL)
    with open(mutfile) as f:
        allMutSmp = f.readline().strip().split("\t")
        resL = [a for a in map(formatSampleName, allMutSmp) if a in resL]
    return resL

# <codecell>

def loadTarIntSmp(gslistfile):
    resD = defaultdict(list)
    with open(gslistfile) as f:
        line = f.readline()
        while line:
            crt_t, crt_s = line.strip().split("\t")
            resD[crt_t] = crt_s.split(";")
            line = f.readline()
    return(resD)

# <codecell>

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

# <codecell>

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

# <codecell>

def loadCNV(cnvfile, smpsL):
    '''
    load cnv matrix; and convert to binary 1/0
    '''
    resD = defaultdict(list)
    with open(cnvfile) as f:
        _, allMutSmp = f.readline().strip().split("\t",1)
        allMutSmp = map(formatSampleName, allMutSmp.split("\t"))
        smpIndx = [id for (id, v) in enumerate(allMutSmp) if v in smpsL]
        resD['gene'] = map(allMutSmp.__getitem__, smpIndx)
        line = f.readline()
        while line:
            crt_g, crt_v = line.strip().split("\t",1)                                  
            temp =  map(lambda x: int(float(x) != 0), crt_v.split("\t"))
            resD[crt_g] = map(temp.__getitem__, smpIndx)
            line = f.readline()
    return resD

# <codecell>

# def sortDictValBySmp(mutD):
#     smp = mutD['gene']
#     del mutD['gene']
#     for k, v in mutD.items():
#         mutD[k] = v.sort(key = dict(zip(v, smp)).get)
#     mutD['gene'] = sorted(smp)
#     return mutD

# <codecell>

def jaccardCoeff(lista, listb):
    cnti = 0.0; cntj = 0.0; cntmatch = 0.0 
    if len(lista) == 0 and len(listb) == 0:
        return 0
    else:
        for i,j in zip(lista, listb):
            if i !=0:
                cnti = cnti + 1
                if j != 0:
                    cntj = cntj + 1
                    if i == j:
                        cntmatch = cntmatch + 1
    if (cnti + cntj) > 0 :    
        return cntmatch / (cnti + cntj)        
    else :
        return 0
    


# <codecell>

tarRegD = loadTarReg(keyRegSumfile)
tarIntSmpD = loadTarIntSmp(gslistfile)
mutExpsmpL = getMutExpSmp(expfile, mutfile)
expD = loadExp(expfile, mutExpsmpL)

# <codecell>

cnvD = loadCNV(cnvfile, mutExpsmpL)
mutD = loadMut(mutfile, mutExpsmpL)
smpCnv = cnvD['gene']
smpMut = mutD['gene']

# <codecell>

mutD['gene'] = mutSmp; cnvD['gene'] = cnvSmp

# <codecell>

def mergeMutCnv(mutD, cnvD):
    '''
        sum up mutation matrix with cnv matrix, binary sum up, 
        1-0, 0-1, 1-1 ==> 1; 0-0 ==> 0
        sample order follo
    '''
    mutSmp = mutD['gene']; cnvSmp = cnvD['gene']
    sample = list(set(mutSmp).intersection(set(cnvSmp)))
    del mutD['gene'];del cnvD['gene']

    
    mutGene = mutD.keys(); cnvGene = cnvD.keys()
    gene = list(set(mutGene).intersection(set(cnvGene)))
    mutOnlyGene = list(set(mutGene) - set(gene))
    cnvOnlyGene = list(set(cnvGene) - set(gene))
    
    mutDF = pd.DataFrame(mutD, index=mutSmp).T
    cnvDF = pd.DataFrame(cnvD, index=cnvSmp).T
    
    #should concatenating two list union of all genes
    vartDF = mutDF.loc[gene,sample] + cnvDF.loc[gene, sample]  
    vartDF = concat([vartDF, mutDF.loc[mutOnlyGene,sample], \
                     cnvDF.loc[cnvOnlyGene, sample]])
    
    
#     print mutDF.loc[gene[100:105],sample[20:25]]
#     print cnvDF.loc[gene[100:105],sample[20:25]]
#     print vartDF.loc[gene[100:105],sample[20:25]]
    
    
    mutD['gene'] = mutSmp; cnvD['gene'] = cnvSmp
    
    vartDF = vartDF.where(vartDF == 0, 1)
    vartD = vartDF.T.to_dict(outtype='list')
    vartD['gene'] = sample
    return vartD

# <codecell>

vartD = mergeMutCnv(mutD, cnvD)

# <codecell>

import matplotlib.pyplot as plt
jacdCoeffD = defaultdict(float)
for k in mutD.keys():
    jacdCoeffD[k] = jaccardCoeff(mutD[k], cnvD[k])

# <codecell>


from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages(figd + '/hist_JaccardCoeff_somaticWcnv_070821014.pdf')

fig = plt.figure()
# plot Jaccard Coefficient
n, bins, pathces = plt.hist(jacdCoeffD.values(), 50, \
                            normed=1,facecolor = 'g', alpha = 0.75)
plt.xlabel('Jaccard Coeff')
plt.title('Histogram of Jaccard coeff \n between SM and CNV for each gene')
# plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
plt.axis([0, 0.55, 0, 90])
plt.ylabel('Percentage %')
plt.grid(True)
# plt.show()

# plt.savefig(pp, format='pdf')
pp.savefig(fig)

## plot CNV samples
fig = plt.figure()
cnvSmp = cnvD['gene']
del cnvD['gene']
n, bins, pathes = plt.hist(map(sum, cnvD.values()), 50,\
                            normed = 1, facecolor = 'r', alpha = 0.75)
plt.xlabel('CNV samples')
plt.ylabel('Percentage %')
plt.title('Hist of CNV sample ')
plt.grid(True)
plt.show()
# plt.savefig(pp, format='pdf')
pp.savefig(fig)

## plot somatic mutation samples
fig = plt.figure()
mutSmp = mutD['gene']
del mutD['gene']
n, bins, pathes = plt.hist(map(sum, mutD.values()), 50,\
                            normed = 1, facecolor = 'b', alpha = 0.75)
plt.xlabel('Somatic Mutation samples')
plt.ylabel('Percentage %')
plt.title('Hist of Somatic Mutation sample ')
plt.grid(True)
plt.show()
# plt.savefig(pp, format='pdf')
pp.savefig(fig)

pp.close()

mutD['gene'] = mutSmp
cnvD['gene'] = cnvSmp

# <codecell>

####----main function debugging mode

# for tgene in tarRegD.keys()[:1]: 
tgene = "GRK5"
# def startOpt4Gene(tgene, tarIntSmpD, tarRegD, mutD, expD, outputDir):
tIntSmp = tarIntSmpD[tgene]
allRegsL = tarRegD[tgene]
intMutSmpIdL = [id for (id, s) in enumerate(vartD['gene']) if s in tIntSmp]
intExpSmpIdL = [id for (id, s) in enumerate(expD['gene']) if s in tIntSmp]

## the regMutD here includes both mut and cnv
regMutD = {k:map(v.__getitem__, intMutSmpIdL) for (k,v) in vartD.items() if k in allRegsL}

regExpD = {k:map(v.__getitem__, intExpSmpIdL) for (k,v) in expD.items() if k in allRegsL}

expMutRegL = set(regExpD.keys()).intersection(set(regMutD.keys()))
outTempMut = output + "_" + tgene + "_regMut"
outTempExp = output + "_" + tgene + "_exp"

outTempMutH = open(outTempMut,'w')
outTempExpH = open(outTempExp,'w')

outTempMutH.write('gene\t'+"\t".join( map(mutD['gene'].__getitem__, intMutSmpIdL)) + "\n")
for k,v in regMutD.items():
    if k in expMutRegL: 
        outTempMutH.write(k + "\t" + "\t".join(map(str,v)) + "\n")

debugSmp = expD['gene']
outTempExpH.write('gene\t')
for i in intExpSmpIdL:
#     print debugSmp[i]
    outTempExpH.write( debugSmp[i] + "\t" )
outTempExpH.write("\n")

# outTempExpH.write('gene\t'+"\t".join(map(expD['gene'].__getitem__, intExpSmpIdL)) + "\n")

## target expession in the first row

debugSmp = expD[tgene]
outTempExpH.write(tgene + '\t')
for i in intExpSmpIdL:
    print i
    print debugSmp[i]
#     outTempExpH.write( debugSmp[i] + "\t" )
outTempExpH.write("\n")

# outTempExpH.write(tgene + "\t" + "\t".join(map(str,map(expD[tgene].__getitem__, intExpSmpIdL))) +"\n")

for k,v in regExpD.items():
    if k in expMutRegL:
        outTempExpH.write(k + "\t" + "\t".join(map(str,v)) + "\n")  
outTempMutH.close()
outTempExpH.close()

outTemp = output + "/" + tgene + ".tsv"
#     cmd = "/ifs/home/c2b2/ac_lab/jh3283/tools/R/R-3-02/bin/Rscript \
#     /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step2-2_regKeyRegulators_v5.r\
#     --vanilla --exp " +  outTempExp + " --mut " + outTempMut + " --output " + outTemp 


# cmd = "~/tools/R/R_current/bin/Rscript /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step3-4_greedyOptCorr.r \
# --vanilla --exp " +  outTempExp + " --mut " + outTempMut + " --output " + outTemp 
# print cmd
# subprocess.Popen(cmd, shell = True)



#     rtncode  = subprocess.call(cmd, shell = True)
#     if rtncode != 0 :
#         print "Error in regression!"
#         sys.exit()
#     print tgene + " Done"

# <codecell>

print len(set(mutSmp).intersection(set(cnvSmp)))
print len(set(mutSmp))
print len(set(cnvSmp))
          

# <codecell>


# <codecell>


# <codecell>


# <codecell>

    

# <codecell>

### test module 2 
xx = [2,0,3]
xxlab = ['b','a','c']
yy = [0,2,6]
xD = {'xx' : xx, 'yy': yy}

from pandas import DataFrame
import pandas as pd
print df
df = pd.DataFrame(xD, index=xxlab).T
xx = [0,1,0,0]
xxlab = ['b','a','d','c']
yy = [0,0,1,0]
xD = {'xx' : xx, 'yy': yy}

df2 = pd.DataFrame(xD, index=xxlab).T
# print df2
mycol = ['a', 'b']
# print df[mycol]
# print df2[mycol]
df3 = df.loc[['xx','yy'],mycol] + df2.loc[['xx','yy'], mycol]
print df3
# print df.T.to_dict(outtype='list')
df3 = df3.where(df3 == 0, 1)
# print df3
print df3.loc[['xx', 'yy'],['a','b']]
df3['zz'] = -1
concat([df3, df3.loc[['xx'],:]])

# <codecell>



# <codecell>




# <codecell>



# <codecell>

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

# <codecell>


# <codecell>


# <codecell>


# <codecell>


# <codecell>


# <codecell>


# <codecell>


# <codecell>


# <codecell>



# <codecell>


# <codecell>


# <codecell>

def format2Dlist(l2d):
    out = []
    for i in l2d:
        out.append("["+",".join(map(str,i)) + "]")
    return ";".join(out)

# <codecell>



# <codecell>


# <codecell>


# <codecell>


# <codecell>


