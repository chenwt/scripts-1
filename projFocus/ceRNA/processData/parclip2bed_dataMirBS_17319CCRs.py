#!/usr/bin/python
#J.HE
'''
input: PAR-CLIP microRNA binding site
output : bed file like binding site
purpose: data process, final
last modified: Jun 27, 2014
'''

import os, sys, re
from collections import defaultdict

argv = sys.argv[1:]

input = argv[0]
output = argv[1]


outDict = defaultdict(list)
with(open(input)) as f:
    line = f.readline()
    while line:
        line = re.split(r"\||:|-|\[", \
                        line.strip().replace("]","").replace("chr","")) 
        key = "-".join(line[1:])
        outDict[key].append(line[0]) 
        line = f.readline()

fhout = open(output, 'w')
for k,v in sorted(outDict.items()):
    if len(v) > 1:
        outRec = k.replace("-", "\t") + "\t" + ",".join(v)
    else:
        outRec = k.replace("-", "\t") + "\t" + v[0]
    fhout.write(outRec + "\n" )  

fhout.close()
print "[END]"
