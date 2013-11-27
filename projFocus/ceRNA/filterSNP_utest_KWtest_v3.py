#!/usr/bin/python
# -*- coding: utf-8 -*-
#J.HE
#input:  -s <snp annotation file> 
#        -e <exp annotation file> 
#        -o <output file name> 
#        -j <number:cutoff for pval adjust 1e-6>
#output:    <test statistic file: Gene-SNP pair with p-value for U/KW test>
#           <gene-snp pair snp GT file>
#           <log file >


import os
import sys, getopt
import math
import re
from scipy import stats
import numpy as np

#define a gene class
class Gene:
  def __init__(self, name, chrom, pos, strand, val):
    self.name = name
    self.chrom = chrom
    self.pos = pos
    self.strand = strand
    self.val = val
    
  def __lt__(self, other):
    return self.pos < other.pos

argv = sys.argv[1:]
inps      = ''
inpe      = ''
fdr_cut	  = 1e-6 
outp      = ''
try:
   opts,args = getopt.getopt(argv,"hs:e:j:o:",["snpfile=","expfile","fdradj=","ofile="])
except getopt.GetoptError:
   print 'python filterSNP_utest_KWtest_v3.py -s <snp annotation file> -e <exp annotation file> -o <output file name> -j <number:cutoff for pval adjust 1e-6>' 
   sys.exit(2)

for opt, arg in opts:
  if opt == '-h':
   	print 'python filterSNP_utest_KWtest_v3.py -s <snp annotation file> -e <exp annotation file> -o <output file name> -j <number:cutoff for pval adjust 1e-6>'
   	sys.exit()
  elif opt in ("-s","--snpfile"):
     inps = arg
  elif opt in ("-e","--expfile"):
     inpe = arg
  elif opt in ("-j","--adjpcutoff"):
     fdr_cut = arg
  elif opt in ("-o","--ofile"):
     outp = arg
if inps == '' or inpe == '' or outp == '':
    print "input parameter ERROR!"
    print 'python filterSNP_utest_KWtest_v2.py -s <snp annotation file> -e <exp annotation file> -o <output file name> -j <number:cutoff for pval adjust 1e-6>'
    sys.exit(2)

region_cut	  = 1000000
chroms	  = range(1,23) + ['X','Y','x','y'] + map(str,range(1,23))
fdr_cut = float(fdr_cut)
outlog = outp + ".log"
outlogf	  = open(outlog,'w')


print('Input file:\t' + inps)
print('Database file:\t' + inpe)
print('Output file:\t'+ outp + ".adjPass_" + str(fdr_cut) + ".mat" + "\t" + outp + ".adj.snp")
print('Log file:\t'+ outlog)

outlogf.write('Input file:\t'    + inps	  + "\n")
outlogf.write('Database file:\t' + inpe	  + "\n")
outlogf.write('Output file:\t'   + outp	  + "\n")
outlogf.write('Log file:\t'      + outlog + "\n")
outlogf.flush()
##----------------------------
#useful functions
def chrom2Num(chrom):
  if chrom == 'X' or chrom == 'x':
    chrom = 23
  if chrom == 'Y' or chrom == 'y':
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

def isInRange(x,y,cut):
  return math.fabs(int(x)-int(y)) < int(cut)

def rUTest(x,y):
  if len(x) != len(y):
    print "ERROR: different lenght for Utest"
    sys.exit(2)
  else:
    p = r.wilcox_test(x,y)['p.value']
    z = r.wilcox_test(x,y)['statistic']['W']
  return([z,'%e' % float(p)])
def kwTest(xlist,ylist):
  x = np.array(xlist)
  y = np.array(ylist)
  [t,p] = stats.mstats.kruskalwallis(x,y)
  return([t,'%e' % float(p)])

def pAdj(plist,type):
  plist = map(float,plist)
  plist_sorted = sorted(plist)
  m = len(plist)
  p_new = []
  if   type == 'b':
    for i, pval in enumerate(plist):
        p_new.append(min(pval * m / (plist_sorted.index(pval) + 1),1))
  elif type == 'f':
    """ need to develop """
    pass
  return(p_new) 

def binarySearch(chromPosArray, pos):
  low = 0 - 1
  high = len(chromPosArray) - 1
  while high >= low:
    mid = (low + high) / 2
    gene = chromPosArray[mid]
    if gene.pos > pos:
      high = mid - 1
    elif gene.pos < pos:
      low = mid + 1
    else:
      return mid + 1  
    #if gene.pos == pos:
    #    return mid
    #elif gene.pos < pos:
    #    low = mid + 1
    #else:
    #    high = mid - 1
  return low

# load all gene's tss information with expression data
###----loading_Gene_Information----
print "loading gene expression data..."
allChroms = map(str,range(24)) + ['X','Y','x','y']
chromArray = [[] for _ in range(24)] ##sort chrom pos information
genePosDict  = {} ## stor chrom_pos genename, values mapping information
i = 0
with open(inpe, 'r') as inpef:
  line = inpef.readline()  
  headerExp = line.strip().split("\t",4)[4]
  line = inpef.readline()
  while len(line) > 0:
    i = i + 1
    [identifier,chrom,pos,strand,val] = line.strip().split("\t")[:5]
    line = inpef.readline()
    if chrom  in allChroms:
      chrom = chrom2Num(chrom)
      pos = int(pos)
      gene = Gene(identifier, chrom, pos, strand, val)
      chromArray[chrom].append(gene)

for chromPosArray in chromArray:
  chromPosArray.sort()

print "Total Number of Genes input:\t" + str(i) 
outlogf.write( "Total Number of Genes input:\t" + str(i) + "\n")
outlogf.flush()

#####-----------------------loading_snp_data
#SNP_A-1780977   22      30598097        -       2       2       2       2       2
cnt_out = 0 
chromDict = [{} for _ in range(24)]
pvalArray = []

outpf            = open(outp,'w')
snpdataDict      = {}
cnt_out          = 0
snpNameSet = set()

with open(inps, 'r') as inpsf:
  line = inpsf.readline()
  headerSnp = line.strip().split("\t",4)[4]
  line = inpsf.readline()
  while len(line) > 0:
    [snpname, chrom, pos, strand, val ] = line.strip().split("\t",4)
    snpNameSet.add(snpname)
    line = inpsf.readline()

snpdataDict = dict()
snpdataDict.fromkeys(snpNameSet) 

cnt_snp = 0 
cnt_p = 0
with open(inps, 'r') as inpsf:
  line = inpsf.readline()
  line = inpsf.readline()
  while len(line) > 0:
    [snpname, chrom, pos, strand, val ] = line.strip().split("\t",4)
    line = inpsf.readline()
    snpdata = Gene(snpname, chrom, pos, strand, val)
    snpdataDict[snpname] = snpdata
    cnt_snp = cnt_snp + 1 
    if chrom in allChroms: 
      chrom = chrom2Num(chrom)
      tempStart = int(pos) - region_cut
      tempEnd = int(pos) + region_cut
      ## find all genes that have this snp around
      index = binarySearch(chromArray[chrom], tempStart)
      if index < len(chromArray[chrom]) :
        expgene = chromArray[chrom][index]
        while tempEnd > expgene.pos:
          cnt_out = cnt_out + 1
          snpval = map(int, val.split("\t"))
          expval = map(float,expgene.val.split("\t"))
          genename = expgene.name
          #[zscore,pval] = rUTest(exp,snpval)
          [zscore,pval] = kwTest(expval,snpval)
          cnt_p = cnt_p + 1
          if cnt_p % 1000 == 0 :
              print "testing snp_gene pair:\t" + str(cnt_p) + "..."
          pvalArray.append([genename,snpname,zscore,pval])
          index = index + 1
          if index < len(chromArray[chrom]):
            expgene = chromArray[chrom][index]
          else:
            break
print ("Total Number of SNP input: " + str(cnt_snp))
print ("Total Number of pval: " + str(cnt_p))
outlogf.write("Total gene_snp pair tested:\t" + str(cnt_out) + "\n")
outlogf.flush()
#------------p_value_correction_and_output----
print "doing pvalue correction...."
outpf	= open(outp + ".adjPass_" + str(fdr_cut) + ".mat",'w')
outpsnpf= open(outp + ".adj.snp",'w')

p_old	= [row[3] for row in pvalArray]
p_new	= pAdj(p_old,'b')
[outdata[0].append(outdata[1]) for outdata in zip(pvalArray,p_new) if outdata[1] < fdr_cut]
outpf.write("Gene" + "\t" + "snpname" + "\t" + headerSnp + "\n")
outpsnpf.write("Gene" + "\t" + "snpname" + "\t" + "statistic" + "\t" + "pval" + "\t" + "pval_adj" + "\n")

cnt_out = 0 
print "writing output...."
for out in pvalArray:
    outpsnpf.write("\t".join(map(str,out)) + "\n")
    if len(out) > 4 :
        cnt_out = cnt_out + 1
        if cnt_out % 1000 == 0:
            print "writing output gene snp pair \t" + str(cnt_out) +"..."
        outpf.write(out[0] + "\t" + out[1] + "\t" + snpdataDict[out[1]].val + "\n")
    else:
        continue
outpsnpf.close()
outpf.close()
print("Total gene_snp pair passed fdr " + str(fdr_cut) + "\t" + str(cnt_out) )
outlogf.write("Total gene_snp pair passed fdr " + str(fdr_cut) + "\t" + str(cnt_out) + "\n")
print "#-------------END------"