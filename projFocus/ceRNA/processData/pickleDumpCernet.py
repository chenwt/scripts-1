#!/usr/bin/python
#J.HE
#Desp.: save data as python binary file, serilizing python processing.

usage = 'python " + sys.argv[0] + " -i <input>  -o <output>'
example = 'python " + sys.argv[0] + " -i <input>  -o <output>'

import sys, os, getopt
import pickle 
import generalUtils as gu
input = sys.argv[1:][0]
print input
output = input + ".pickle"

cnt = 1
with open(input) as f:
    line = f.readline()
    while line:
        if cnt == 1:
            net = gu.Cernet(line.strip().split("\t",2)[:2])
        else:
            net.__append__(line.strip().split("\t",2)[:2])
        cnt = cnt + 1
        line = f.readline()
outputH = open(output, 'w')

pickle.dump(net, outputH)

print "#----END------"
