#!/usr/bin/env python
#input: 1.<file:expression data of genes>
#	2.<file: snp GT file>
#	3.<file: meth mval file>
#	4.<file: the full path of cnv level 3 files> 

import os
import sys,getopt


#arg = sys.argv[1:]

inpe = 'test_exp_chr22.mat'
inps = 'test_snp_chr22.mat'
inpm = 'test_meth_chr22.mat'

###----mainDict to store the location of points
mainDict = {} 
with opne(inpe) as inpef:
  for i,line in enumerate(inpef):
    if i==0 :
      pass
    else:
      [id,chr,pos,strand,val] = line.strip().split("\t")   
      key = chr + ":" + pos  
      mainDict[key] = []

print "Number of Genes: " + str(len(mainDict.keys()))
###-----dict to store snp points values
with open(inps) as inpsf:
  for i, line in enumerate(inpsf):
    if i == 0:
      pass
    else:
      [id,chr,pos,strand,val] = line.strip().split("\t")   
      key = chr + ":" + pos  
###-----dict to store meth points values
with open(inps) as inpsf:
  for i, line in enumerate(inpsf):
    if i == 0:
      pass
    else:
      [id,chr,pos,strand,val] = line.strip().split("\t")   
      key = chr + ":" + pos  
###-----get Gene CNV infor 
## need more time for this part since 


#---setting up parameters
amp_cut = 0.3
del_cut = -0.3
#----for_homo_cnv---
#amp_cut = 0.7
#del_cut = -0.9
#----end------------

