#!/usr/bin/python
#J.HE
#Desp: to convert TCGA somatic mutation Level2 .maf file into a matrix of gene * sample mutate/not
#input: <1. file: TCGA somatic mutation level2. maf> 
#       <2. file: sample name want to include in final output file>
#output: <1: file: somatic mutation indicator file of gene by sample>
#       <2: stats about the the somatic mutation distribution>
 
USAGE = "python selectTargetRegionRow.py \n\
        -i \n\
        -t \n\
        -c \n\
        -o "
EXAMPLE =""

import os,sys
import getopt
import re
from sets import Set

class Region():
    def __init__(self, chr, ps, pe):
        self.chr = chr
        self.ps = int(ps)
        self.pe = int(pe)
    def overlap(self, reg2):
        if self.chr == reg2.chr:
            if self.pe < reg2.ps or self.ps > reg2.pe:
                return 0
            else:
                return 1
        else:
            return 0

def gsCont(gene,sample):
    return(gene + ":" + sample)

def chr2Num(x):
    if x in ['x','X']:
        x = 23
    elif x in ['Y','y']:
        x = 24
    elif x not in range(23) + ['x','X','Y','y']:
        x = 25
    return(int(x)) 

argv = sys.argv[1:]
try:
    opts,args = getopt.getopt(argv,"hi:t:c:o:")
except getopt.GetoptError:
    print "-i <input maf file> -o <output matrix file> "
    sys.exit()

for opt, arg in opts:
    if opt == "-h":
        print ""
    elif opt == "-i":
        mafFile = arg
    elif opt == "-t":
        targetfile = arg
    elif opt == "-c":
        REGION_CUT = int(arg)
    elif opt == "-o":
        outFile = arg

#load target region
tReg = [[] for _ in range(26)] 
cnt = 0  
with open(targetfile) as f:
    line = f.readline()
    while line:
        if re.match("^#", line) or re.match("^gene",line):
            cnt = cnt + 1
            pass
        else:
            genename, chrom, pstart, pend, strand  =\
                    re.split("\t", line.strip())[:5]
            tempReg = Region(chrom, pstart, pend)
            tempReg.strand = strand
            tempReg.genename = genename
            tReg[chr2Num(chrom)].append(tempReg)
            cnt = cnt + 1
        line = f.readline()
print "total target region\t" + str(cnt) 

for tRegEle in tReg:
    tRegEle.sort(key=lambda x:x.ps)

#load mutation genes and samples
mutInfo = Set()
cntLine = 0 
allGene = Set()
sampleName = Set()
outFileHander = open(outFile,'w')
cnt = 0
cntline = 0
with open(mafFile) as f:
    line = f.readline()
    while line:
        if re.match("^[#|gene]",line) or cnt ==0:
            outRecord = line
            outFileHander.write("geneTarge\tchrTarge\tpsTarget\tpeTarget\tgeneMut\tpsMut\tpeMut\t" +
                                line.split("\t",1)[1])
            cnt = cnt +1
        else:
            tempgene, tempchr, tempps, tempend, tempval = \
                    line.strip().split("\t",4)
            temprs = long(tempps) - REGION_CUT
            tempre = long(tempend) + REGION_CUT
            mutReg = Region(tempchr, temprs, tempre)
            mutReg.genename = tempgene
            for reg in tReg[chr2Num(tempchr)]:
                if reg.overlap(mutReg):
                    outRecord = [reg.genename, reg.chr, reg.ps, reg.pe, \
                    tempgene, tempps, tempend, tempval] 
                    cnt = cnt +1
                    outFileHander.write("\t".join(map(str,outRecord)) + "\n")
                else:
                    pass 
        cntline = cntline + 1
        line = f.readline()
print "lines wrote\t" + str(cnt)
print "lines process\t" + str(cntline)

outFileHander.close()
print "#----END----"
                
