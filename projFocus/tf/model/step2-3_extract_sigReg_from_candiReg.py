#!/usr/bin/python
#J.HE
'''
input: result file of candiReg which output all regulators
output: input.sigreg
'''
import os, sys

input = sys.argv[1:][0]
pvalCut = 0.1
out = []; cnt = 0
with(open(input)) as f:
    line = f.readline()
    while line:
        if line[0] == "#" :
            out.append(line)
        else:
            g, coeff, pvalue = line.strip().split("\t")
            if (float(coeff) != 0 and float(pvalue) <= pvalCut) \
               and g[0] != "(":
                out.append(line)
            else:
                cnt = cnt + 1
        line = f.readline()

if cnt >0 :
    info1, val, = out[2].strip().split("\t")
    out[2] = "\t".join([info1, str(int(val) - cnt)])  + "\n" 

fout = open(input + ".sigreg", 'w')
for i in out:
    fout.write(i)
fout.close()
# print "[---END----]"
