#!/usr/bin/env python
# -*- coding: utf-8 -*-
#input: 1. <file:file names of the level 3 cnv file> 2. <file: genes,and location of interest:format: identifier: chr. pos. strand>
#output: 1. <file: .mat.anno file including all cnv info for the input files>
#Desp.: given genename, output CNV value in each gene
#Sample  Chromosome      Start   End     Num_Probes      Segment_Mean
#BITES_p_TCGAb61_SNP_S_GenomeWideSNP_6_C01_697154        1       61735   82650   11      0.5307
      

usage = "python step1-3_getGeneCNVMatLevel2.py \
        -f <file:full path of cnv level3 files,one each line> \
        -g <file: gene file with position inform; gene,chr,start,end,strand> \
        -o <filename: output file name>"

import os
import sys, getopt

argv = sys.argv[1:]


try:
  opts,args = getopt.getopt(argv, "h:f:g:o:")
except:
  print usage  
  sys.exit()
for opt,arg in opts:
  if opt == '-h':
    print usage  
    sys.exit()
  elif opt == '-f':
    inpc = arg
  elif opt == '-g':
    inpg = arg
  elif opt == '-o':
    out = arg
    outlog = out + ".log"

##-----setting parameters
amp_cut	   = float(0.2) ## 1.22 copy
del_cut	   = float(-0.2) ## 0.819 copy
nesp	   = 1  ##number of lines to escape at the begining of each cnv file
nval	   = 5 - 1 ###which column has the value, 5 for TCGA level3 data

class Gene:
  def __init__ (self, chr, start, end, strand):
    self.chr   = chr
    self.start = int(start)
    self.end   = int(end)
    self.strand = strand
    self.cnv = 0.0
  def updateCNV(self,cnval):
     cnval = float(cnval)
     if abs(cnval) > abs(self.cnv):
         self.cnv = cnval
     else:
         self.cnv 
  def initCNV(self):
      self.cnv = int(0)
  def notOverlap(self, region):
    ''' compare whethe first gene is in the second '''
    if self.chr == region.chr :
	if self.start > region.end or self.end < region.start:
	   res = True
	else:
	   res = False
    else:
      res = True
    return res
  def printAll(self):
      attrs = vars(self)
      print ', '.join("%s: %s" % item for item in attrs.items())

class CNV:
    def __init__ (self, chr, start, end, val):
        self.chr   = chr
        self.start = int(start)
        self.end   = int(end)
        self.strand = '+'
        val = float(val)
        if val < amp_cut and val > del_cut:
            self.val = 0.0
        else :
            self.val = val 
    def printAll(self):
        attrs = vars(self)
        print ', '.join("%s: %s" % item for item in attrs.items())
    
##--------load all cnv filenames 
outlogf = open(outlog,'w')
fnArray = []
with open(inpc) as inpf:
  for line in inpf.readlines():
    fnArray.append(line.strip())
    
nCNVSamples = len(fnArray)
print "CNVsample number \t" + str(nCNVSamples)

###----loading_Gene_Information----
geneList = []
with open(inpg) as inpgf:
  for line in inpgf.readlines():
      [identifer,chrom,start,end,strand] = line.strip().split("\t")
      gene = Gene(chrom,start,end,strand)
      gene.name = identifer
      geneList.append(gene)
nGene = len(geneList)
outlogf.write("number_of_genes\t" + str(nGene) + "\n")
print "number of genes:\t" +str(nGene)
#####--------------------------------loading CNV data. 
idfile = -1
outputH = open(out,'w')
for fn in fnArray:
  try:
    fntempf = open(fn)
    idfile = idfile + 1
    print fn
    if idfile ==0:
        outputH.write("barcode"+"\t" + '\t'.join([gene.name for gene in geneList]) + "\n")             
    for i, line in enumerate(fntempf):
      if i < nesp:
        pass
      else:
        [identifier,chrom,pStart,pEnd,probeNum,val] = line.strip().split("\t")
        cnv = CNV(chrom,pStart,pEnd,val) 
        cnv.sample = fn
        tempcnt = 0 
        for gene in geneList:
            if  gene.notOverlap(cnv) == True :    
                continue
            else:  
                gene.updateCNV(cnv.val)
    fntempf.close()  
    outputH.write(fn+"\t" + '\t'.join(map(str,[gene.cnv for gene in geneList])) + "\n")
    for gene in geneList: gene.initCNV()              
  except IOError as e:
    print "I/O error({0}): {1}".format(e.errno, e.strerror)
    sys.exit()
    
outlogf.write("number_of_cnv_level3_files\t" + str(idfile + 1) + "\n")  

outputH.close()
outlogf.close()
print "#--------------DONE--------"
