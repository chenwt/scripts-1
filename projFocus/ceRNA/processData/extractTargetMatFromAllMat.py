#!/usr/bin/env python
#J.HE
#Desp: for projFocus ceRNA, given caner target genes, extract ceRNET regulator genes from giving ceRNET network
#Input: -i -d -o
usage = "~/tools/python/Python_current/python \
/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/extractTargetMatFromAllMat.py \
-i <genelist> \
-d <big matrix > \
-o <output>" 
example = 'python \
/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/step1-1_extractCenetRegulator.py \
 -i /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/gslist/CG_combine_02172014.list \
 -d /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/cnv/brca_cnv_l3_tumor_Mar-23-2014.matrix \
 -o /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/gslist/target_cnv.matrix'

##----------------------------
#functions
def getElement(a,l):
	[res for ele in l if ele==a]
	return res
##----------------------------

import os,sys,getopt
import re

argv = sys.argv[1:]
inputfile  = ''
inputDB = ''
output   = ''
tyep = '' 

try:
    opts,args = getopt.getopt(argv,"hi:d:t:o:")
except getopt.GetoptError:
    print usage + "\n" + example
    sys.exit()

for opt, arg in opts:
    if opt == '-h':
        print usage +"\n" + example
        sys.exit()
    elif opt in ("-i"):
        inputfile = arg
    elif opt in ("-d"):
        inputDB = arg
    elif opt in ("-t"):
        type = arg
    elif opt in ("-o"):
        output = arg

print inputfile
print inputDB
print output

target = []
with open(inputfile) as f:
    for line in f.readlines():
	    tempGene = line.strip()
	    target.append(tempGene)

print "Input gene number:\t" + str(len(target))
cnt = 0 
outputH = open(output,'w')
with open(inputDB) as f:
    line = f.readline()
    while line:
        gene,val=line.split("\t", 1) 
        if re.match("^barcode|gene", gene) or re.findall("\.01A",line):
            outputH.write(line)
        elif gene in target:
            outputH.write(gene + "\t" + val)
            if type == 'uniq':
                target.remove(gene)
            cnt = cnt + 1
        else:
            pass
        line = f.readline()
outputH.close()

print "Output records number:\t" + str(cnt)
print "#---END---"

