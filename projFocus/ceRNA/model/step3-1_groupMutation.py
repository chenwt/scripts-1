#!/usr/bin/python
#J.HE
#Desp.: given the mutation matrix, group mutations based on gene/GeneRegion


import sys,getopt
from mutationRecord import MutRecord
from collections import defaultdict 
import numpy as np
import pickle 
argv = sys.argv[1:]
input = ''
output = ''
usage = 'python ' + sys.argv[0] + ' -i <input>  -o <output>'
example = 'python ' + sys.argv[0] + ' -i <input>  -o <output>'
try:
    opts,args = getopt.getopt(argv,"hi:o:")
except getopt.GetoptError:
    print usage + "\n" + example 
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print usage + "\n" + example 
        sys.exit()
    elif opt in ("-i"):
        input = arg
    # elif opt in ("-o"):
        # output = arg
outputFreq = input  + ".groupByGene.freq"
outputMat  = input + ".groupByGene.binary.matrix"
outputMapfile = input + ".groupByGene.groupMap.pickle"

print('Script path:\t'+ sys.argv[0])
print('Input file:\t' + input)
print('Output file:\t'+ output)


mutDict = defaultdict(list)
mutValDict = defaultdict(list)
with(open(input)) as f:
    head = f.readline()
    line = f.readline()
    while line:
        tgene, chrom, ps, pe, val =\
                line.strip().split("\t", 4)
        crtMut = MutRecord(tgene, tgene, chrom, ps, pe, val)
        if mutDict.get(tgene,''):
            mutDict[tgene].append(mutDict)
            mutValDict[tgene] = mutValDict[tgene]  + \
                  np.array(val.split("\t"),dtype=int ) 
        else:
            mutDict[tgene].append(mutDict)
            mutValDict[tgene] = np.array( val.split("\t"), dtype=int)
        line = f.readline()

print len(mutValDict.items())
outputFreqH = open(outputFreq, 'wt')
outputMatH  = open(outputMat, 'wt')

outputFreqH.write("geneTarget\trecurrence\n")
outputMatH.write('Gene\t'+ "\t".join(head.strip().split("\t",4)[4:]) + "\n" )
for k,v in mutValDict.items():
    outputFreqH.write(k+"\t"+str(sum(v))+"\n")
    outputMatH.write(k + "\t" + "\t".join(map(lambda x:str(min(x,1)), v)) + "\n") 
pickle.dump([mutDict,mutValDict], open(outputMapfile,'w'))
