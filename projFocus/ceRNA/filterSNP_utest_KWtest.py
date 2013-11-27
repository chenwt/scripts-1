#!/usr/bin/python
#J.HE
#input: <file:snp.annot file> <file: exp annot file>
#output: <file: Gene-SNP pair with p-value for U/KW test>

import os
import sys, getopt
import math
import re
from rpy import *
from scipy import stats
import numpy as np

argv = sys.argv[1:]
inps = ''
inpe = ''
outp = ''
try:
  opts,args = getopt.getopt(argv,"hs:e:o:",["snpfile=","expfile'","ofile="])
except getopt.GetoptError:
  print 'python filterSNP_utest_KWtest.py -s <snp annotation file> -e <exp annotation file> -o <output file name>' 
  sys.exit(2)

for opt, arg in opts:
  if opt == '-h':
  	print 'python linkSNO_EXO.py -s <snp annotation file> -e <exp annotation file> -o <output file name>'
  	sys.exit()
  elif opt in ("-s","--snpfile"):
    inps = arg
  elif opt in ("-e","--expfile"):
    inpe = arg
  elif opt in ("-o","--ofile"):
    outp = arg
if inps == '' or inpe == '' or outp == '':
	print "input parameter ERROR!"
	print "python linkSNO_EXO.py -s <snp annotation file> -e <exp annotation file> -o <output file name> "
	sys.exit(2)

cutoff	  = 1,000,000
chroms	  = range(1,23) + ['X','Y','x','y'] + map(str,range(1,23))
fdr_cut	  = 1e-8  
outlogf	  = open(outlog,'w')
outp = outp + "_" + str(fdr_cut)
outlog = outp + ".log"

print('Input file:' + inps)
print('Database file:' + inpe)
print('Output file:'+ outp)
print('Log file:'+ outlog)

outlogf.write('Input file:\t'    + inps	  + "\n")
outlogf.write('Database file:\t' + inpe	  + "\n")
outlogf.write('Output file:\t'   + outp	  + "\n")
outlogf.write('Log file:\t'      + outlog + "\n")


# load all gene's tss information with expression data
print "loading gene expression data..."
chrArrayExp = [{} for _ in range(24)]
with open(inpe) as inpef:
	for i, line in enumerate(inpef): 
		if i == 0:
			continue
		else:
			[genename,chrom,pos,strand,val] = line.strip().split("\t",4)
			if chrom in chroms :
				tempP = pos.split(":")
				if len(tempP) > 1 :
					tempP = map(int,tempP)
					if  strand == "+" :
						pos = min(tempP)
					elif strand == "-" :
						pos = max(tempP)
					else:
						print "ERROR! Strand information error!"
						break
				else :
					pos = int(pos)

				if chrom == "X" :
					chrom = 23
				if chrom == "Y":
					chrom = 24
				chrom = int(chrom)
				if chrom in range(1,25) :
					chrArrayExp[chrom - 1][pos]=[genename,strand,val]
				else :
					continue
			else:
				continue
outlogf.write( "Total Number of Genes input:\t" + str(i+1) + "\n")

#----------useful_functions-----------
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
  m	= len(plist)
  if   type == 'b':
    p_new = map(lambda x:x * m, plist)
  elif type == 'f':
    """ need to develop """
    pass
  return(p_new) 

#--------------------load snp data 
print "loading snp data..."
outpf = open(outp,'w')
outpf.write("Gene\tSNP\tZscore\tpval\n")
pvalArray = []
snpdataDict = {}
cnt_out	    = 0 
with open(inps) as inpsf:
	for i, line in enumerate(inpsf):
		if i == 0:
			header_snpSample = line.strip().split("\t",4)[3]
			continue
		else:
			[snpname,chrom,pos,strand,val] = line.strip().split("\t",4)
			snpdataDict[snpname] = [chrom,pos,strand,val]		
			if chrom in chroms:
				pos = int(pos)
				if chrom == "X" :
					chrom = 23
				if chrom == "Y":
					chrom = 24
				chrom = int(chrom)
				if chrom in range(1,25) :
					for key in chrArrayExp[chrom - 1].keys():
						if isInRange(pos,key,cutoff):
							cnt_out  = cnt_out + 1
#							print "pos : \t" + str(pos) + "\t in range of \t" + str(key)
							exp_data      = chrArrayExp[chrom - 1][key][2]
							snpval	      = map(int,val.strip().split("\t"))
							exp	      = map(float,re.compile("\s+").split(exp_data.strip()))
							genename      = chrArrayExp[chrom - 1][key][0]
							#[zscore,pval] = rUTest(exp,snpval)
							[zscore,pval] = kwTest(exp,snpval)
							#print "zscore:" + str(zscore) + "\t" + "pvalue" + str(pval) 
							pvalArray.append([genename,snpname,zscore,pval])
							outpf.write(chrArrayExp[chrom - 1][key][0]  \
								+ "\t" + snpname + "\t" \
								+ str(zscore) + "\t" + str(pval) + "\t" \
								+ "\n")
				else :
					continue
			else:
				continue

outpf.close()
#----------p_value_correction_and_output----
outpf = open(outp".adjPass_" + str(fdr_cutt),'w')
p_old	= [row[3] for row in pvalArray]
p_new	= pAdj(p_old,'b')
[outdata[0].append(outdata[1]) for outdata in zip(pvalArray,p_new) if outdata[1] < fdr_cut]

outpf.write("Gene" + "\t" + "snpname" + "\t" + header_snpSample + "\n")
for out in outdata:
  if len(out) > 4 :
    outpf.write(outdata[0] + "\t" + outdata[1] + "\t" + "\t".join(snpdataDict[outdata[1]]) + "\n")
  else:
    continue

outpf.close()
#os.system('~/tools/R/R_current/bin/Rscript /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/pvalAdj.r ' + outp  )
#os.system('cp ' + outp + ' ' + outp +'.bak'  )
#print "running r adjustment"
#outpf = open(outp,'w')
#with open(outp + ".adjusted") as inpe:
#  for i,line in enumerate(inpe):
#    if i ==0:
#      pass
#    else:
#      [gene,snpname,val] = line.strip().split("\t",2)
#      outpf.write(snpname,"\t", "\t".join(snpdataDict[snpname]) + "\n")
#print str(i+1) + "SNPs kept"
#outpf.close()
print "#--------------DONE-------"
