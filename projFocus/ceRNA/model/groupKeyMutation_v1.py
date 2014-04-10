#!/usr/bin/python
#J.HE

import re
import sys, os 
from collections import defaultdict
from mutationRecord import MutRecord  

# CWD='/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/'

keyRegfile = "/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/candiReg/temp-gslist.Gintset_Mar31.txt_42/CDH11_candidateRegs_Mar-31-2014.txt"
mutPromoterfile = "~/DATA/projFocus/data/03102014/tcgal2som/genome.wustl.edu__Illumina_All.maf.matrix.promoter2k.Mar-20-2014.matrix.sorted"
mutMirBSfile = "/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/sigMut/brca_som_cupidBS_Apr-7-2014_v1.matrix"
gene = 'CDH11'
debug = 'False'
smpSizeCut = 3

keyRegInfoDict = defaultdict()
regList = [] 
with open(keyRegfile) as f:
    line = f.readline()
    while line :
        if re.match("^#",line) and not re.findall("^regulator", line):
            info, val = line.strip().replace("#","").split("\t",1)
            keyRegInfoDict[info]= val  
        else:
            reg, coeff, pval = line.strip().split("\t",2)
            regList.append(reg)
            assert float(coeff) > 0. 
        line = f.readline()

r2PvalCut = 0.05
if int(keyRegInfoDict['sigReg']) < 1:
    sys.exit()
if float(keyRegInfoDict['r2.pval']) > r2PvalCut:
    sys.exit()

if debug:
    print keyRegInfoDict['sigReg']
    print regList


def extractMut(regs, mutfile):
    '''
    extract regulators mutation information from big mutation matrix
    more details
    '''
    mutList = []
    with(open(mutfile)) as f:
        head = f.readline()
        line = f.readline()
        while line:
            targ, mutg, tarchr, mutps, mutpe, val = \
                    line.strip().split("\t",5)
            if mutg in regList or targ in regList:
                crtMut = MutRecord(targ, tarchr, mutg, mutps, mutpe, val) 
                crtMut.freq = val.strip().split("\t")
                mutList[mutg].append(crtMut)
            line = f.readline()
    return mutList

mutAllList  = extractMut(regList, mutMirBSfile) + extractMut(regList, mutPromoterfile)

## get gint sample

## check recurrenct
# get all mutations by occurence 
outMut = []

if crtMut.freq >= smpSizeCut:
    outMut.append(crtMut)
    del mutAllDict[crtMut]
## grouping 
setCov = 0 
setCovOverlap = 0
setSize = 1  
'''
setCov >=0, <= len(Gint sample) numbe of smps have mutations in one of the gene
setCovOverlap, #of smps have mutations in more than a gene in of set
setSize = 1, #of genes in the set  <= len(regList)
find a set M of k genes(Argmax setSize) with maximum setCov
'''


