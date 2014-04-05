#!/usr/bin/python
#J.HE
#Desp.: transform fasta format to id:seq format


import sys, os, re

fasta = sys.argv[1:][0]
out   = sys.argv[1:][1]

outH = open(out, 'w')

cnt = 0
with(open(fasta)) as f:
    line = f.readline()
    while line:
        if re.match("^>", line):
            pt = r"\s+|="
            g, _, rv = re.split(pt, line.strip())[:3]
            g = g.split("_",2)[2]
            coord = rv.replace("chr","") 
            if cnt != 0 :
                outH.write("\n")
            outH.write(g +"\t" + coord + "\t")
        else:
            outH.write(line.strip())
        cnt = cnt + 1
        line = f.readline()

print "%s lines process!" % cnt

