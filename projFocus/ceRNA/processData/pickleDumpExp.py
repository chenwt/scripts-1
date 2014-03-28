#!/usr/bin/python
#J.HE
#Desp.: load expression data, and save data into a gene_by_chr dict format, 
# as python binary file, serilizing python processing.

usage = 'python " + sys.argv[0] + " -i <input>  -o <output>'
example = 'python " + sys.argv[0] + " -i <input>  -o <output>'

import sys, os, getopt
import pickle 
from operator import itemgetter
from collections import defaultdict, OrderedDict

input = sys.argv[1:][0]
annof = sys.argv[1:][1]
print input
output = input + ".pickle"
outputlog = output + ".log" 
outputlogH = open(outputlog, 'w')
##-----------------load gene annotation infor
# annoDict = defaultdict(list)
# with open(annof) as f:
#     line = f.readline()
#     while line:
#         g, chrom, tss, tse, strand = line.strip().split("\t")
#         annoDict[g] = [chrom, tss, tse, strand]
#         line = f.readline()
# print len(annoDict)
# outputH = open(annof + ".pickle", 'w') 
# pickle.dump(annoDict, outputH)

annoDict = pickle.load(open(annof + ".pickle")) 

##-----------------load expression
geneByChr = defaultdict(dict) 
cnt = 1
with open(input) as f:
    line = f.readline()
    while line:
        if cnt == 1:
            allSamples = line.strip().split("\t")
            geneByChr['samples']={'all': allSamples}
        else:
            tmpGene,vals = line.strip().split("\t",1)
            try :
                geneByChr[annoDict[tmpGene][0]].update(\
                        {tmpGene:vals.split("\t") })
            except IndexError:
                outputlogH.write(tmpGene + "\n")
        cnt = cnt + 1
        line = f.readline()
print geneByChr.keys()
outputH = open(output, 'w')
pickle.dump(geneByChr, outputH)
outputH.close()
outputlogH.close()

print "#----END------"
