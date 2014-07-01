#!/usr/bin/python
#J.HE
'''
usage: $scriptname <inputfilename> <outputfilename>
input: downloaded data of mirBinding site interaction with Argo
output: bed file similar with each binding site as one row
purpose: data process- final
lastModidfied: Jun 27 2014

'''

import os, re, sys
from collections import defaultdict

argv = sys.argv[1:]

input = argv[0]
output = argv[1]

outDict = defaultdict(list)
with(open(input)) as f:
    line = f.readline()
    if re.findall(r"chrom", line, re.I):
        line = f.readline()
    while line:
        line = line.strip().replace("chr","")
        chr, ps, pe, _, score, strand, mir, _ = line.split("\t", 7)
        key = '-'.join([chr,ps,pe,strand])
        outDict[key].extend([mir, float(score)]) 

        if len(outDict[key]) < 3:
            line = f.readline()
            continue

        outDict[key][-1] = max(outDict[key][-3], outDict[key][-1])
        del outDict[key][-3]

        line = f.readline()

print "uniq BS", len(outDict.keys())

fhout = open(output, 'w')
for k, v in outDict.items():
    outRec = k.replace("-","\t") + "\t" +  \
             str(v[-1]) + "\t"+ ",".join(v[:-1])  
    fhout.write( outRec + "\n")

print "[END]" 
