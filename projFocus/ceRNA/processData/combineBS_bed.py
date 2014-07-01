#!/usr/bin/python
#J.HE
'''
input: 4 bed files of mirBinding site
output : 1 big bed file with mir BS
purpose: data prcess; complete

'''
import os, sys, re, time
from collections import defaultdict

tstart = time.time()
argv = sys.argv[1:]

fclash = argv[0]
fparclip = argv[1]
fmc2004 = argv[2]
fcupid = argv[3]

output = argv[-1]

outRecDict = defaultdict(list) 

with(open(fclash)) as f:
    line = f.readline()
    while line:
        line = line.strip()
        chr, ps, pe, _ = line.strip().split("\t",3) 
        key = chr + "\t" + ps + "\t" + pe
        outRecDict[key].append("clash") 
        line = f.readline()

print fclash + " loaded" 

with(open(fparclip)) as f:
    line = f.readline()
    while line:
        line = line.strip()
        chr, ps, pe, _ = line.strip().split("\t",3) 
 
        key = chr + "\t" + ps + "\t" + pe
        outRecDict[key].append("parclip") 

        line = f.readline()

print fparclip + " loaded" 

with(open(fmc2004)) as f:
    line = f.readline()
    while line:
        line = line.strip()
        chr, ps, pe, _ = line.strip().split("\t",3) 
        
        key = chr + "\t" + ps + "\t" + pe
        outRecDict[key].append("mc2004") 

        line = f.readline()
print fmc2004 + " loaded" 

fhout = open(output, 'w')

cnt = 0 
with(open(fcupid)) as f:
    line = f.readline()
    while line:
        line = line.strip()
        chr, ps, pe, _ = line.strip().split("\t",3) 
        
        key = chr + "\t" + ps + "\t" + pe
        if not outRecDict.get(key,0):
            outRec = key + "\t" + "cupid;"          
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
print "total uniq bind site", cnt 
print "time elapsed: :", time.time() - tstart, "s" 
print "[END]"
