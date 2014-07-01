#!/usr/bin/python
#J.HE
'''
input: tf binding and mir binding 
output : 1 big bed file with all regulatory BS
purpose: data prcess; complete

'''
import os, sys, re, time
from collections import defaultdict

tstart = time.time()
argv = sys.argv[1:]

finput1 = argv[0]
finput2 = argv[1]

output = argv[-1]

outRecDict = defaultdict(list) 

with(open(finput1)) as f:
    line = f.readline()
    while line:
        line = line.strip()
        chr, ps, pe, info = line.strip().split("\t",3) 
        key = chr + "\t" + ps + "\t" + pe
        outRecDict[key].append(info) 
        line = f.readline()

print finput1 + " loaded" 

fhout = open(output, 'w')

cnt = 0 
with(open(finput2)) as f:
    line = f.readline()
    while line:
        line = line.strip()
        chr, ps, pe, info = line.strip().split("\t",3) 
        
        key = chr + "\t" + ps + "\t" + pe
        if not outRecDict.get(key,0):
            outRec = key + "\t" + info   
        else :
            outRec = key + "\t" + ";".join(outRecDict[key])
            del outRecDict[key]
        fhout.write(outRec  + "\n")
        cnt = cnt + 1
        line = f.readline()

for key,v in outRecDict.items():
    outRec = key + "\t" + ";".join(list(set(v)) ) 
    fhout.write(outRec  + "\n")
    cnt = cnt + 1

fhout.close()

print "total uniq bind site", cnt 
print "time elapsed: :", time.time() - tstart, "s" 
print "[END]"
