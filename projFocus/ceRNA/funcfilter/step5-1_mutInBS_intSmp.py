#!/usr/bin/python
#J.HE
'''
input: 1. mut matrix, with header including the sample name
       2. gslist, target gene and it's intact samples
       3. tar-reglist, result from either lasso or greedy, target and regulator
 output: bed file, 5 column, last colums is the mutated in total intact samples,
and the total number of regulated target genes

Comments: 
'''

import sys,getopt
from collections import defaultdict
import re
from generalUtils import formatBarcode

argv = sys.argv[1:]
input = ''
output = ''
usage = 'python ' + sys.argv[0] + ' -i <bed file of mutation list> \
-o <output> -s <gene intact sample list> -g <target gene ceRNA driver list> \
        -m <mutation matrix> -o <output bed with 6 column> '
example = 'python ' + sys.argv[0] + ' -i <input>  -o <output>'
try:
    opts,args = getopt.getopt(argv,"hs:g:m:o:")
except getopt.GetoptError:
    print usage + "\n" + example 
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print usage + "\n" + example 
        sys.exit()
    elif opt in ("-s"):
        fintsmp = arg
    elif opt in ("-m"):
        fmut = arg
    elif opt in ("-g"):
        ftarreg = arg
    elif opt in ("-o"):
        output = arg
# print('Script path:\t'+ sys.argv[0])
# print('Input file:\t' + input)
# print('Output file:\t'+ output)


rTarDict = defaultdict(list)

with(open(ftarreg)) as f:
    line = f.readline()
    while line:
        line = re.split("\t|;", line.strip())
        for r in line[1:]:
            rTarDict[r].append(line[0]) 
        line = f.readline()

print "target regulator file loaded"

intsmpDict = defaultdict(list)
with(open(fintsmp)) as f:
    line = f.readline()
    while line:
        line = re.split("\t|;", line.strip())
        intsmpDict[line[0]] = line[1:]
        line = f.readline()

print "intact sample file loaded"

outHeader = "\t".join(["genename","chrom", "mutps", "mutpe", \
                       "mutsmp", "totalTarget"])

fhout = open(output,'w')
fhout.write(outHeader + "\n")

def getUnionOfIntSmp(regulator, rTarDict, tarIntSmpDict):
    intSmpUnion = []
    for t in rTarDict[regulator]:
        intSmpUnion.extend(tarIntSmpDict[t])
        # print regulator + "\t" + t + "\t" + str(len(tarIntSmpDict[t]))
    return list(set(intSmpUnion))

geneArry = []
with(open(fmut)) as f:
    line = f.readline()
    allsmpArry = map(formatBarcode, line.strip().split("\t") )
    line = f.readline()
    while line:
        line = line.strip().split("\t")
        crtIntSmp = getUnionOfIntSmp(line[5], rTarDict, intsmpDict)
        idx = [i for (i,s) in enumerate(allsmpArry) if s in crtIntSmp] 

        if len(idx) == 0:
            line = f.readline()
            continue

        mutSmpCnt = reduce(lambda x,y: x+y, map(int,map(line.__getitem__, idx)))
        totalTarCnt = len(rTarDict[line[1]])

        outRec = line[1:3] + line[6:8] + [str(mutSmpCnt), str(totalTarCnt)]
        fhout.write("\t".join(outRec) + "\n") 

        line = f.readline()
fhout.close()

print "#[---END---]"
