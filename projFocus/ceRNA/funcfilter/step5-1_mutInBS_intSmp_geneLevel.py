#!/usr/bin/python
#J.HE
'''
input: 1. mut matrix, with header including the sample name
       2. gslist, target gene and it's intact samples
       3. tar-reglist, result from either lasso or greedy, target and regulator
 output: bed file, 5 column, last colums is the mutated in total intact samples,
and the total number of regulated target genes

Comments:  compare to the first version, collect all int sample for one single
gene before and then count the mutation number

Status: underwork
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

def getUnionOfIntSmp(regulator, rTarDict, tarIntSmpDict):
    intSmpUnion = []
    for t in rTarDict[regulator]:
        intSmpUnion.extend(tarIntSmpDict[t])
    return list(set(intSmpUnion))

geneArry = []
regIntSmpDict = defaultdict(list)
regMutDict = defaultdict(list)

with(open(fmut)) as f:
    line = f.readline()
    allsmpArry = map(formatBarcode, line.strip().split("\t") )
    line = f.readline()
    while line:
        line = line.strip().split("\t")
        gname = line[5]
        crtIntSmp = getUnionOfIntSmp(gname, rTarDict, intsmpDict)
        regIntSmpDict[gname].extend(crtIntSmp)
        regIntSmpDict[gname] = list(set(regIntSmpDict[gname]))
        line = f.readline()

outRecDict = defaultdict(list)
regMutSmpCnt = defaultdict(list)
with(open(fmut)) as f:
    line = f.readline()
    line = f.readline()
    while line:
        line = line.strip().split("\t")
        gname = line[5]
        crtIntSmp= regIntSmpDict[gname]

        idx = [i for (i,s) in enumerate(allsmpArry) if s in crtIntSmp] 
        if len(idx) == 0:
            line = f.readline()
            continue

        mutSmpCnt = reduce(lambda x,y: x+y, map(int,map(line.__getitem__, idx)))
        regMutSmpCnt[gname].append(mutSmpCnt)
        totalTarCnt = len(rTarDict[gname])

        outRecDict[gname] = [gname,str(totalTarCnt)]
        line = f.readline()


#output
outHeader = "\t".join(["genename","totalTarget", "mutCnt", "mutsmp" ])

fhout = open(output,'w')
fhout.write(outHeader + "\n")


for k,v in outRecDict.items():
        crtCnt = regMutSmpCnt[k]
        if len(crtCnt) > 1:
            outRec = v + [str(len(crtCnt)), str(reduce(lambda x,y:x+y, crtCnt))]
        else:
            outRec = v + ['1'] + [str(crtCnt[0])]
        fhout.write("\t".join(outRec) + "\n") 

fhout.close()
print "#[---END---]"
