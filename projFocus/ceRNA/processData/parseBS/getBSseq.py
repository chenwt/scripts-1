#!/usr/bin/python
#J.HE
#Desp.: given the cupid predict, and sequence, return cupid predict sequence

import os, sys
import re
from collections import defaultdict

import sys,getopt
argv = sys.argv[1:]
input = ''
output = ''
usage = 'python ' + sys.argv[0] + ' -i <input>  -o <output>'
example = 'python ' + sys.argv[0] + ' -i <input>  -o <output>'
try:
    opts,args = getopt.getopt(argv,"hp:s:o:")
except getopt.GetoptError:
    print usage + "\n" + example 
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print usage + "\n" + example 
        sys.exit()
    elif opt in ("-p"):
        cupidfile = arg
    elif opt in ("-s"):
        seqfile = arg
    elif opt in ("-o"):
        output = arg
print('Script path:\t'+ sys.argv[0])
print('Input file:\t' +cupidfile)
print('Output file:\t'+ output)

seq3pUTR = defaultdict(list)
cnt = 0
with(open(seqfile)) as f:
    line = f.readline()
    while line:
        if not re.match("^Symbol", line):
            _, genename, refseqid, _, seq  = re.split(r":|\s+",line.strip())
            seq3pUTR[genename].append(seq)
            cnt = cnt+1
        line = f.readline()
print " %s sequence for %s genes" % (cnt, len(seq3pUTR.keys()) )

geneprev = ''
outseq = ['Symbol', 'sequence']
with(open(cupidfile)) as f:
    line = f.readline()
    while line:
        if not re.match("^Symbol", line):
            _, genename, refseqid, _, ps, pe = \
                    re.split(r":|\s+|,", line.strip().replace("[","").replace("]",""))
            # if not geneprev == genename:
                # outseq = list(set(outseq))
                # print geneprev + "\t" + ";".join(outseq)
                geneprev = genename
                outseq = []
                for gSeq in  seq3pUTR[genename] :
                    ps = max( int(ps) - 1 - 25, 0)
                    pe = min( int(pe) + 25, len(gSeq) )  
                    if not gSeq[ps:pe] in outseq:
                        outseq.append(gSeq[ps:pe]) 
            else:     
                for gSeq in  seq3pUTR[genename] :
                   ps = max( int(ps) - 1 - 25, 0)
                   pe = min( int(pe) + 25, len(gSeq) )  
                   if not gSeq[ps:pe] in outseq:
                       outseq.append(gSeq[ps:pe]) 
                geneprev = genename
        line = f.readline()

print geneprev + "\t" + ";".join(outseq)
