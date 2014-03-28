#!/usr/bin/python
#J.HE
#Desp.: save data as python binary file, serilizing python processing.

usage = 'python " + sys.argv[0] + " -i <input>  -o <output>'
example = 'python " + sys.argv[0] + " -i <input>  -o <output>'

import sys, os, getopt
import pickle 
input = sys.argv[1:][0]
print input
output = input + ".pickle"

data = []
with open(input) as f:
    line = f.readline()
    while line:
        data.append(line.strip().split("\t"))
        line = f.readline()

for line in data:
    print line
outputH = open(output, 'w')
pickle.dump(data, outputH)

print "#----END------"
