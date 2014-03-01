#!/usr/bin/python
#J.HE
#Desp: to convert somatic mutation Level2 .maf file into a matrix of mutation * sample maf
#input: <1. file: somatic mutation .maf file generated using step-6_get_MutMaf.py > 
#       <2. file: sample name want to include in final output file>
#output: <1: file: somatic mutation map file of gene by sample>
#       <2: stats about the the somatic mutation distribution>

 
import os,sys
import getopt
from sets import Set

argv = sys.argv[1:]
try:
    opts,args = getopt.getopt(argv,"hm:o:")
except getopt.GetoptError:
    print ""
    sys.exit(2)

for opt, arg in opts:
    if opt == "-h":
       print ""
    elif opt == "-m":
        mafFile = arg
    elif opt == "-o":
        outFile = arg

def gsCont(gene,sample):
     return(gene + ":" + sample)

#load mutation genes and samples
mutInfo = {}
cntLine = 0 
allMut = Set()
sampleName = []
with open(mafFile) as f:
    f.readline()
    for line in f.readlines():
        lineTemp = line.strip().split("\t")
        [mut, sample] = [gsCont(lineTemp[0],lineTemp[1]), lineTemp[3]]
	mutInfo.update({gsCont(mut,sample):float(lineTemp[2])})
        allMut.add(mut)
	if lineTemp[3] not in sampleName:
	  sampleName.append(lineTemp[3])
	else :
	  pass
        cntLine = cntLine + 1

print "nonZero:\t" + str(len(mutInfo))
print "total Mut nunber:\t" +  str(len(allMut))
print "Total sample number:" + str(len(sampleName))
outFileHander = open(outFile,'w')
outFileHander.write('gene' + "\t" + "\t".join(sampleName) +"\n")

cntMut = 0
for mut in allMut:
    cntMut = cntMut + 1
    if cntMut % 1000  == 0:
        print "prcessed mutation\t" + str(cntGene)  
    geneMutTemp = [] * len(sampleName)
    for sample in sampleName:
	temp = gsCont(mut,sample)
        if temp in mutInfo.keys():            
            geneMutTemp.append(mutInfo[temp])
        else:
            geneMutTemp.append(0)
    outFileHander.write(mut + "\t" + '\t'.join(map(str,geneMutTemp)) + "\n")

outFileHander.close()
sys.exit()
print "#----END----"
                
        
    
              
            
    
         
