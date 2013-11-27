#!/usr/bin/env python
#J.HE
#input:	<file: list of mutations>
#outpu: <file: matrix of mutaions by patients>

import os
import sys
import operator
import pandas as pd


inp  = 'ALL.NoAlt_0_fisher_exact.txt'
outp = 'ALL.NoAlt_0_fisher_exact.txt.matrix'
itemIndexList = [0,1,2,8,10,21]

dataArray = {}
with open(inp) as inpf:
  for i, line in enumerate(inpf):
    if i == 0:
      header = line.strip().split()
    else:
      [pidType, chrom, pos, gene, func, pvalue] = \
        operator.itemgetter(*itemIndexList)(line.strip().split())
      pidType = pidType[:10]
      print [gene,pidType,func,pvalue]
      #sys.exit()
      dataArray[chrom + "_" + pos]= [gene,pidType,func,pvalue]
      if i == 5:
          break


dataDF = Series(dataArray,index= ['chrom_pos', 'gene', 'pidType', 'func', 'pvalue'])
dataDF.groupby("pidType")['pidType'].transform('mean')
dataNew = dataDF['pvalue']
dataNew.index =  pd.MultiIndex.from_arrays([df['chrom_pos'], dataDF['pidType'])
print dataNew.index