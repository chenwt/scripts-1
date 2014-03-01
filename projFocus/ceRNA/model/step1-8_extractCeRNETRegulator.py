#!/usr/bin/env python
#J.HE
#Desp: for projFocus ceRNA, given caner target genes, extract ceRNET regulator genes from giving ceRNET network
#Input: -i -d -o

usage = '~/tools/python/Python_current/python step1-8_extractCeRNETRegulator.py  -i <input target gene>  -d <cancer specific CeRNET_from Hua-sheng> -o <output file name> '
example = 'python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/validSNP.py -i -d  -o '

##----------------------------
#functions
def getElement(a,l):
	[res for ele in l if ele==a]
	return res
##----------------------------

import os,sys,getopt
import re

argv = sys.argv[1:]
inputfile  = ''
inputDB = ''
output   = ''


try:
  opts,args = getopt.getopt(argv,"hi:d:o:")
except getopt.GetoptError:
  print usage + "\n" + example
  sys.exit()

for opt, arg in opts:
    if opt == '-h':
          print usage +"\n" + example
          sys.exit()
    elif opt in ("-i"):
          inputfile = arg
    elif opt in ("-d"):
          inputDB = arg
    elif opt in ("-o"):
          output = arg

print inputfile
print inputDB
print output

target = []
with open(inputfile) as f:
    for line in f.readlines():
	tempGene = line.strip().split()[0] if re.findall(r'\s+',line) else line.strip()
	target.append(tempGene)
print "Input targets:\t " + str(len(target))

targetInDB = []
cnt = 0 
with open(inputDB) as f:
	for line in f.readlines():
	        [rna1, rna2] = line.strip().split("\t")[:2]
		if (rna1 in target) or (rna2 in target):
			targetInDB.extend([rna1, rna2])
			cnt = cnt + 1
		else:
			continue
targetInDB = list(set(targetInDB))
print "Input ceRNET interaction:\t" + str(cnt)
print "target ceRNET regulators:\t" + str(len(targetInDB))

outputH = open(output,'w')
outputH.write("\n".join(targetInDB)+ "\n")
outputH.close()
print "#---END---"

