#!/usr/bin/env python
#J.HE
#Desp: generated to link somatic mutation with diff expressed genes, using somMut region infor
#input: 1. <file:anno somatic mutation .mat.anno file>  2. <file: annot exp.mat file> 
#output: 1. <file: somticbyExp.anno.mat>      

import os
import sys, getopt
import numpy as np
import pickle

argv = sys.argv[1:]
usage = "usage: python mergeExpSomMutation.py  -s <file:full path of cnv level3 files,one each line> -e <file: gene_exp.mat> -o <filename: output file name>"
example = "python mergeExpSomMutation.py -s /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/data/somaticMutation/brca_som_selectedSample_level2.mat.5colanno -e /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/exp/brca_exp_l3_731_DEG.mat.singleTSS.anno -o brca_somForDeg.mat  "
try:
  opts,args = getopt.getopt(argv, "hs:e:o:")
except:
  print usage 
  print example
  sys.exit()

for opt,arg in opts:
  if opt == '-h':
    print usage 
    print example
    sys.exit()
  elif opt == '-s':
    inpSomFile = arg
  elif opt == '-e':
    inpExpFile = arg
  elif opt == '-o':
    outFile = arg
    outlog = outFile + ".log"
outlogf = open(outlog,'w')

##-----setting parameters
region_cut = 1000000 # 1M
#-------end-

def chrom2Num(chrom):
  if chrom == 'X':
    chrom = 23
  if chrom == 'Y':
    chrom = 24
  chrom = int(chrom) - 1
  return(chrom)

def num2Chrom(num):
  num = int(num)
  if num == 22:
     num = 'X' 
  elif num == 23:
    num = 'Y' 
  else:
    num = num + 1
  return(str(num))

###----loading_Gene_Information----
allChroms = map(str,range(22)) + ['X','Y']
chromArray = [[] for _ in range(24)] #keep all position for each chrom
genePosDict  = {} #keep all chr, pos identifier information
with open(inpExpFile) as f:
  for i,line in enumerate(f):
    if i == 0:
      headerExp = line.strip()
    else:
      [identifer,chrom,pos,strand] = line.strip().split("\t")[:4]
      if chrom in allChroms: 
	chrom = chrom2Num(chrom)
	if genePosDict.get((chrom,pos),0): 
	  genePosDict[(chrom,pos)].append([identifer])
	else:
	  genePosDict[(chrom,pos)] = [identifer]
	chromArray[chrom].append(int(pos))
      else:
	continue
outlogf.write("number_of_genes\t" + str(i) + "\n")

print "total number of chr_pos \t" + str(len(genePosDict.keys()))

#####--------------------------------loading somMutation data. 
chromDict = [{} for _ in range(24)]
expSomLink = []
outFileHandler = open(outFile,'w')
with open(inpSomFile) as f:
  for i, line in enumerate(f):
     if i == 0 :
       pass
       print "skip header line: \t" + line[:50] + "\t..."
       outFileHandler.write("Gene\tsomMutation\t"+ line.strip().split("\t",5)[5] + "\n") 
     else:
       [somGene,chrom,start,end,strand,val] = line.strip().split("\t",5)
       if chrom in allChroms:
          chrom = chrom2Num(chrom)
          tempStart = int(start) - region_cut
          tempEnd = int(end) + region_cut
          #---find all TSS site for one somaticMutation gene 
	  for p in chromArray[chrom]:
	    if p > tempStart and p < tempEnd :
	      expGene = genePosDict[(chrom,str(p))]
	      if len(expGene) > 1:
		expSomLink.extend(map(list,zip(expGene,somGene * len(expGene),val * len(expGene))) )
	      else:
		expSomLink.append([expGene,somGene,val])
       else:
	  continue

print "somatic mutation processed:\t" + str(i)
print "output binary linkfile..."

#pickle.dump(expSomLink, open('expSomLink.p','wb'))
#expSomLink = pickle.load(open('expSomLink.p','rb'))

###----- sumup somatic mutation for each expGene
cntGene = 1
for expGene in genePosDict.values():
  temp = [map(int,row[2].split("\t")) for row in expSomLink if row[0] == expGene] 
  if len(temp) > 0 and len(temp[0]) > 0: 
    tempNP = np.array(temp)
    sumTempNP = np.sum(tempNP,0)
    cntGene = cntGene + 1
    outFileHandler.write(str(expGene[0]) +"\tsomMutation\t"+ "\t".join(map(str,sumTempNP.tolist())) + "\n")
  else:
    continue

outFileHandler.close()
print "DEG genes output:\t" + str(cntGene) 
print "#------END-----"
    
