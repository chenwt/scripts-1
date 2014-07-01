#!/usr/bin/python
#J.HE
'''
input: TFBS file from hua-sheng
output: bed file like version
purpose: data process
status: underwork 
'''

import os,sys, re

input, output = sys.argv[1:3]

fhout = open(output, 'w')

with(open(input)) as f:
    line = f.readline()
    if re.findall(r"(?i)target", line):
        line = f.readline()
    while line:
        tf, target, region, score = line.strip().replace("chr", "").split()
        chr, ps, pe = re.split("-|:", region)
        outRec = "\t".join([chr, ps, pe, score + "|" + tf + "|" + target]) 
        fhout.write(outRec + "\n")
        line = f.readline()

fhout.close()
print "[END}"
