#!/usr/bin/python
#J.HE
#input: <file:snp.annot file> <file: exp annot file>
#output: <file: Gene-SNP pair with p-value for U/KW test>

import os
import sys, getopt
import math
import re
from scipy import stats
import numpy as np

argv = sys.argv[1:]
inps      = ''
inpe      = ''
fdr_cut	  = 0.5  
outp      = ''
try:
   opts,args = getopt.getopt(argv,"hs:e:j:o:",["snpfile=","expfile","fdradj=","ofile="])
except getopt.GetoptError:
   print 'python filterSNP_utest_KWtest_v2.py -s <snp annotation file> -e <exp annotation file> -o <output file name> -j <number:cutoff for pval adjust 1e-6>' 
   sys.exit(2)

for opt, arg in opts:
  if opt == '-h':
   	print 'python filterSNP_utest_KWtest_v2.py -s <snp annotation file> -e <exp annotation file> -o <output file name> -j <number:cutoff for pval adjust 1e-6>'
   	sys.exit()
  elif opt in ("-s","--snpfile"):
     inps = arg
  elif opt in ("-e","--expfile"):
     inpe = arg
  elif opt in ("-j","--ofile"):
     fdr_cut = arg
  elif opt in ("-o","--ofile"):
     outp = arg
if inps == '' or inpe == '' or outp == '':
    print "input parameter ERROR!"
    print 'python filterSNP_utest_KWtest_v2.py -s <snp annotation file> -e <exp annotation file> -o <output file name> -j <number:cutoff for pval adjust 1e-6>'
    sys.exit(2)

# region_cut	  = 1000000
region_cut	  = 100000000
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

##----------------------------
#useful functions
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


# load all gene's tss information with expression data
###----loading_Gene_Information----
print "loading gene expression data..."
allChroms = map(str,range(22)) + ['X','Y','x','y']
chromArray = [[] for _ in range(24)] ##sort chrom pos information
genePosDict  = {} ## stor chrom_pos genename, values mapping information
with open(inpe) as inpef:
  for i,line in enumerate(inpef):
    if i == 0:
      headerExp = line.strip().split("\t",4)[4]
    else:
      [identifier,chrom,pos,strand,val] = line.strip().split("\t")[:5]
      if chrom  in allChroms:
        chrom = chrom2Num(chrom)
        pos = int(pos)
        genePosDict[(chrom,pos)]    = [identifier,val]
        chromArray[chrom].append(int(pos))
      else:
        continue
for chromPosArray in chromArray:
  chromPosArray.sort(key=int)

outlogf.write( "Total Number of Genes input:\t" + str(i) + "\n")


#####-----------------------loading_snp_data
#SNP_A-1780977   22      30598097        -       2       2       2       2       2
cnt_out = 0 
chromDict = [{} for _ in range(24)]
pvalArray = []

outpf            = open(outp,'w')
snpdataDict      = {}
cnt_out          = 0
with open(inps) as inpsf:
  for i,line in enumerate(inpsf):
    if i == 0:
      headerSnp = line.strip().split("\t",4)[4]
    else:
      [snpname, chrom, pos, strand, val ] = line.strip().split("\t",4)
      snpdataDict[snpname] = [chrom,pos,strand,val]
      if chrom in allChroms: 
        chrom = chrom2Num(chrom)
        tempStart = int(pos) - region_cut
        tempEnd = int(pos) + region_cut
        ## find all genes that have this snp around
        for genePos in [p for p in chromArray[chrom] if p > tempStart and p < tempEnd]:
                cnt_out = cnt_out + 1
                snpval = map(int, val.split("\t"))
                expval = map(float,genePosDict[(chrom,genePos)][1].split("\t"))
                genename = genePosDict[(chrom,genePos)][0]		#[zscore,pval] = rUTest(exp,snpval)
		[zscore,pval] = kwTest(expval,snpval)
		pvalArray.append([genename,snpname,zscore,pval])
      else:
        continue
outlogf.write("Total gene_snp pair tested:\t" + str(cnt_out) + "\n")

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
        outpf.write(out[0] + "\t" + out[1] + "\t" + snpdataDict[out[1]][3] + "\n")
    else:
        continue
outpsnpf.close()
outpf.close()
outlogf.write("Total gene_snp pair passed fdr " + str(fdr_cut) + "\t" + str(cnt_out) + "\n")
print "#-------------END------"