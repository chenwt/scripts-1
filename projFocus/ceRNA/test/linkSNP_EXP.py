#!/usr/bin/python
#J.HE
#input: <file:snp.annot file> <file: exp annot file>
#output: <file: Gene-SNP pair with p-value for U test>

import os
import sys, getopt
import math
import re
from rpy import *

argv = sys.argv[1:]
inp1 = ''
inp2 = ''
outp = ''
try:
  opts,args = getopt.getopt(argv,"hs:e:o:",["snpfile=","expfile'","ofile="])
except getopt.GetoptError:
  print 'python linkSNO_EXO.py -s <snp annotation file> -e <exp annotation file> -o <output file name>' 
  sys.exit(2)

for opt, arg in opts:
  if opt == '-h':
  	print 'python linkSNO_EXO.py -s <snp annotation file> -e <exp annotation file> -o <output file name>'
  	sys.exit()
  elif opt in ("-s","--snpfile"):
    inp1 = arg
  elif opt in ("-e","--expfile"):
    inp2 = arg
  elif opt in ("-o","--ofile"):
    outp = arg
    outlog = outp + ".log"
if inp1 == '' or inp2 == '' or outp == '':
	print "input parameter ERROR!"
	print "python linkSNO_EXO.py -s <snp annotation file> -e <exp annotation file> -o <output file name> "
	sys.exit(2)

# inp1 = 'test.snp'
# inp2 = 'test.exp'
# outp = 'test_out_link.txt'
# outlog = outp + ".log"

print('Input file:' + inp1)
print('Database file:' + inp2)
print('Output file:'+ outp)
print('Log file:'+ outlog)

cutoff = 1000000
chroms = range(1,23) + ['X','Y']

# load all gene's tss information with expression data
chrArrayExp = [{}] * 24
with open(inp2) as inp2f:
	for i, line in enumerate(inp2f): 
		if i == 0:
			continue
		else:
			[genename,chrom,pos,strand,val] = line.strip().split("\t",4)
			if chrom in chroms :
				chrom = int(chrom)
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
				if chrom in range(1,25) :
					chrArrayExp[chrom - 1][pos]=[genename,strand,val]
				else :
					continue
			else:
				continue
print "Total Number of Genes processed: " + str(i+1)

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

# load 
outpf = open(outp,'w')

outpf.write("Gene\tSNP\tZscore\t-logP\tEffectSize\n")
with open(inp1) as inp1f:
	for i, line in enumerate(inp1f):
		if i == 0:
			continue
		else:
			[snpname,chrom,pos,strand,val] = line.strip().split("\t",4)
			if chrom in chroms:
				chrom = int(chrom)
				pos = int(pos)
				if chrom == "X" :
					chrom = 23
				if chrom == "Y":
					chrom = 24
				if chrom in range(1,25) :
					for key in chrArrayExp[chrom - 1].keys():
						if isInRange(pos,key,cutoff):
							exp_data = chrArrayExp[chrom - 1][key][2]
							snp_data = val
							tempval = map(int,snp_data.strip().split("\t"))
							tempLen = [ i for i in range(len(tempval)) if(tempval[i] >=0)]
							es = len(tempLen)
							g = chrArrayExp[chrom - 1][key][2]
							gene = map(int,re.compile("\s+").split(g.strip()))
							gene_es = [gene[i] for i in tempLen]
							val_es = [tempval[i] for i in tempLen]
							[zscore,pval] = rUTest(gene_es,val_es)
							outpf.write(chrArrayExp[chrom - 1][key][0]  \
								+ "\t" + snpname + "\t" + \
								str(zscore) + "\t" + str(pval) + \
								"\t" + str(es) + "\n")
				else :
					continue
			else:
				continue
outpf.close()

# sysCmd = "sort -k 1 " + outp + " " + "> " + outp+".sorted"
# os.system(sysCmd)
# sysCmd = "rm " + outp 

print "DONE"


