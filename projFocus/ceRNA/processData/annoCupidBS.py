#!/usr/bin/python
#J.HE
#Desp.: load in refseq gene annotation, save as RefGene object for later use
#input: refseq gene
#output: pickle
### UNDER DEVELOPMENT!!!

usage = 'python " + sys.argv[0] + " -i <input>  -o <output>'
example = 'python " + sys.argv[0] + " -i <input>  -o <output>'

import sys, os, getopt
from collections import defaultdict
import pickle   
import re
from  generalUtils import RefGene 
input = sys.argv[1:][0]
annofile = sys.argv[1:][2]
print input
output = input + ".pickle"

data = defaultdict(list)

def line2GeneObj(line):
    refname, tchr, tstrand, ttxs, ttxe, tcdss, tcdse, \
        _, _, _, tgenename, _, _, _, _, tgenesymbol =\
        re.split("\t",line.strip())
    if tstrand == "+":
        tgeneObj = RefGene(tgenename, tchr.replace("chr", ""), \
                      tstrand, ttxs,  tcdse, ttxe, ttxs, tcdss)
        tgeneObj.refseqname = refname
    else:
        tgeneObj = RefGene(tgenename, tchr.replace("chr", ""), \
                      tstrand, ttxs, ttxs, tcdss, tcdse, ttxe)
        tgeneObj.refseqname = refname
    print tgeneObj.name
    return(tgeneObj)  

cnt = 0
with open(annofile) as f:
    line = f.readline()
    while line:
        if not re.match(r"^#", line):
            tempgObj = line2GeneObj(line) 
            data[tempgObj.name].append(tempgObj)
            cnt = cnt + 1
        line = f.readline()

with(open(input)) as f:
    line = f.readline()
    while line:
        ptn = r":|\s|\t|\[|,|\]"
         _, tg, _, _, trbs, trbe, _ = re.split(ptn, line.strip()) 
         if data.get(tg,''):
             objlist =  
    line = f.readline()

print "total number of line %s" % cnt 
for k, v in data.items():
    data[k] = v.sort(key=lambda g:g.tss)

outputH = open(output, 'w')
pickle.dump(data, outputH)

print "#----END------"
