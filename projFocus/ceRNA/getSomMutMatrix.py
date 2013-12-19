#!/usr/bin/python
#J.HE
#Desp: to convert TCGA somatic mutation Level2 .maf file into a matrix of gene * sample mutate/not
#input: <1. file: TCGA somatic mutation level2. maf> 
#       <2. file: sample name want to include in final output file>
#output: <1: file: somatic mutation indicator file of gene by sample>
#       <2: stats about the the somatic mutation distribution>

 
import os,sys
import getopt
from sets import Set

argv = sys.argv[1:]
try:
    opts,args = getopt.getopt(argv,"hm:s:o:")
except getopt.GetoptError:
    print ""
    sys.exit(2)

for opt, arg in opts:
    if opt == "-h":
       print ""
    elif opt == "-m":
        mafFile = arg
    elif opt == "-s":
        sampleFile = arg
    elif opt == "-o":
        outFile = arg

def gsCont(gene,sample):
     return(gene + ":" + sample)

#load mutation genes and samples
mutInfo = Set()
cntLine = 0 
allGene = Set()
with open(mafFile) as f:
   #skip header line
    f.readline()
    for line in f.readlines():
    #for i, line in enumerate(f):            
        lineTemp = line.strip().split("\t")
        [gene, sample] = [lineTemp[0], lineTemp[15][:19]]
        #print gene + "\t" + sample
        #mutInfo.append([gene,sample])
        mutInfo.add(gsCont(gene,sample))
        allGene.add(gene)
        cntLine = cntLine + 1
        #if cntLine ==10:
        #    break
        #print mutInfo
        #sys.exit()

print "nonZero:\t" + str(len(mutInfo))
print "total Gene nunber:\t" +  str(len(allGene))

#get ordered sample names
sampleName = []
with open(sampleFile) as f:
    for line in f.readlines():
        sampleName.append(line.strip())
        
print "Total sample number:" + str(len(sampleName))
outFileHander = open(outFile,'w')
outFileHander.write('gene' + "\t" + "\t".join(sampleName) +"\n")

cntGene = 0
for gene in allGene:
    cntGene = cntGene + 1
    if cntGene % 1000  == 0:
        print "prcessed gene\t" + str(cntGene)  
    geneMutTemp = [] * len(sampleName)
    for sample in sampleName:
        #if [gene,sample] in mutInfo:
        if gsCont(gene,sample) in mutInfo:            
            geneMutTemp.append('1')
            #print sampleName.index(sample)
        else:
            geneMutTemp.append('0')
    #print gene + "\t" + '\t'.join(geneMutTemp)
    outFileHander.write( gene + "\t" + '\t'.join(geneMutTemp) + "\n")

outFileHander.close()
print "#----END----"
                
        
    
              
            
    
         