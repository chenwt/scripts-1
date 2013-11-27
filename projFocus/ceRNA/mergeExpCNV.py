#!/usr/bin/env python
#input: 1. <file:file names of the level 3 cnv file> 2. <file: genes,and location of interest:format: identifier: chr. pos. strand>
#output: 1. <file: .mat.anno file including all cnv info for the input files>
#Sample  Chromosome      Start   End     Num_Probes      Segment_Mean
#BITES_p_TCGAb61_SNP_S_GenomeWideSNP_6_C01_697154        1       61735   82650   11      0.5307
      

import os
import sys, getopt
#import numpy as np
#import collections

argv = sys.argv[1:]
try:
  opts,args = getopt.getopt(argv, "h:c:e:o:")
except:
  print "ERROR:usage: python mergeExpCNV.py  -c <file:full path of cnv level3 files,one each line> -e <file: gene_exp.mat> -o <filename: output file name>"
  sys.exit()
for opt,arg in opts:
  if opt == '-h':
    print "Usage: python mergeExpCNV.py  -c <file:full path of cnv level3 files,one each line> -e <file: gene_exp.mat> -o <filename: output file name>"
    sys.exit()
  elif opt == '-c':
    inpc = arg
  elif opt == '-e':
    inpg = arg
  elif opt == '-o':
    out = arg
    outlog = out + ".log"

outlogf = open(outlog,'w')

##-----setting parameters
region_cut = 1000000 # 1M
amp_cut	   = float(0.3)
del_cut	   = float(-0.3)
nesp	   = 1
nval	   = 5 - 1
#------homo_cnv--
#amp_cut = 0.7
#del_cut = -0.9
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

##load all cnv filenames 
fnArray = []
with open(inpc) as inpf:
  for line in inpf.readlines():
    fnArray.append(line.strip())

nCNVSamples = len(fnArray)
###----loading_Gene_Information----
allChroms = map(str,range(22)) + ['X','Y']
chromArray = [[] for _ in range(24)]
genePosDict  = {}
with open(inpg) as inpgf:
  for i,line in enumerate(inpgf):
    if i == 0:
      headerExp = line.strip()
    else:
      [identifer,chrom,pos,strand] = line.strip().split("\t")[:4]
      if chrom  in allChroms: 
	genePosDict[(chrom,pos)] = identifer
	chrom = chrom2Num(chrom)
	chromArray[chrom].append(int(pos))
      else:
	continue
inpgf.close()
outlogf.write("number_of_genes\t" + str(i) + "\n")

for chromPosArray in chromArray:
  chromPosArray.sort(key=int)
# chromPosArray.sort(key=lambda tup:tup[0])
## print "chrom\t" + str(chromArray.index(chromPosArray)) +"\tgenes:\t" + str(len(chromPosArray)) 


#####--------------------------------loading CNV data. 
###CAUTIONS: MAKE SURE THE INPUT CNV FILE NAMES ARE THE SAME ORDER AS IN THE SNP AND METH FILE
idfile = -1
chromDict = [{} for _ in range(24)]
for fn in fnArray:
  try:
    fntempf = open(fn)
    idfile = idfile + 1
    for i, line in enumerate(fntempf):
      if i < nesp:
        pass
      else:
        [identifier,chrom,pStart,pEnd,probeNum,val] = line.strip().split("\t")
        if float(val) > amp_cut or float(val) < del_cut: #filter out cnv value
	  if chrom in allChroms:
	     chrom = chrom2Num(chrom)
             tempStart = int(pStart) - region_cut
             tempEnd = int(pEnd) + region_cut
             for pts in [p for p in chromArray[chrom] if p > tempStart and p < tempEnd]:
	       if pts in chromDict[chrom].keys():
                   chromDict[chrom][pts].append([pStart,pEnd,val,idfile])
               else:
                   chromDict[chrom][pts] = [[pStart,pEnd,val,idfile]]
	  else:
	    continue
        else:
          continue
    print str(idfile + 1) + "\tcnv level3 file processed..."
    fntempf.close()
  except IOError as e:
    print "I/O error({0}): {1}".format(e.errno, e.strerror)
    sys.exit()

outlogf.write("number_of_cnv_level3_files\t" + str(idfile + 1) + "\n")  

###-------------- output the result
#loop the chromDictArray for each chrom
cnt_out = 0 
def myCol(matrix,i):
  return([row[i] for row in matrix])

outf = open(out,'w')
outf.write("entrezID\tGenomicVariation\t" + "\t".join(fnArray ) + "\n")
print "writing gene_cnv points...."
####TODO: how to deal with one sample have multiple cnv near one gene
for idx,oneChromDict in zip(range(len(chromDict)),chromDict):
  chromID = num2Chrom(idx) 
  #iterate all pos in each gene
  for genePos in oneChromDict.keys():
    for points in oneChromDict[genePos]:
      tempCNVall = [0] * nCNVSamples
      idx_cnv = myCol(oneChromDict[genePos],-1)
      for id_cnv, val_cnv in zip(idx_cnv,myCol(oneChromDict[genePos],-2)):
        tempCNVall[id_cnv] = val_cnv
      cnt_out = cnt_out + 1
    #print str(chromID) +"_" + str(genePos) + "\t" + "cnv"+"\t"+"\t".join(map(str,tempCNVall))
    #outf.write(str(chromID) +"_" + str(genePos) + "\t" + "cnv"+"\t"+"\t".join(map(str,tempCNVall)) + "\n")
    outf.write(genePosDict[(chromID,str(genePos))] + "\t" + "cnv" + "\t" + "\t".join(map(str,tempCNVall)) + "\n")
print "cnvs_gene:\t" + str(cnt_out) 
outlogf.write("cnvs_gene:\t" + str(cnt_out) + "\n")  

outf.close()
print "#--------------DONE--------"
