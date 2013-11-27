#!/usr/bin/python
# input: <file: filtered somatic mutation files>
# output: <file: somatic mutations with fisher exact test pval> 
#J.HE

import os
import math
import linecache
from scipy import stats

inp = "ALL.NoAlt_0.txt"
outp = "ALL.NoAlt_0_fisher_exact.txt"

def fisherExactTest(freqlist):
  case    = map(int,freqlist[0:2])
  normal  = map(int,freqlist[2:] )
  oddsratio, pvalue = stats.fisher_exact([case, normal])
  return [oddsratio,pvalue]

fout = open(outp,'w')
#the header
with open(inp) as inpf:
  for i, line in enumerate(inpf):
        if i == 0:
            fout.write(line.strip() + \
            "\t" + "fisher_pvalue" + "\n")
        else: 
            [caseRef,caseAlt,normRef,normAlt] = line.strip().split("\t")[10:14]
            #print [caseRef,caseAlt,normRef,normAlt]
            [oddratio, pvalue] = fisherExactTest([caseRef,caseAlt,normRef,normAlt])
            #print oddratio,pvalue,
            fout.write(line.strip() + "\t" + str(pvalue) + "\n")    

inpf.close()
fout.close()
