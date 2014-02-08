#!/usr/bin/env python
# -*- coding: utf-8 -*-
#input: 1. <file:file names of the level 3 meth file> 2. <file: genes,and location of interest:format: identifier: chr. pos. strand>
#output: 1. <file: .mat.anno file including all meth info for the input files>
#Sample  Chromosome      Start   End     Num_Probes      Segment_Mean
#BITES_p_TCGAb61_SNP_S_GenomeWideSNP_6_C01_697154        1       61735   82650   11      0.5307
      

import os
import sys, getopt

argv = sys.argv[1:]
usage = "python" + sys.argv[1] + " -f <file:full path of methylation level3 files,one each line> -g <file: gene file with positive inform; gene,chr,start,end,strand> -o <filename: output file name>"
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

# #-----------test----
# inpg = 'cancer.genelist.startend_022014.tsv'
# inpc = 'input.files'
# out  = 'test.out' 
# outlog = out + ".log"
# outlogf = open(outlog,'w')
# ##-------test--end----

##-----setting parameters
nesp	   = 2  ##number of lines to escape at the begining of each meth file
###-----functions------
def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False
    
class Gene:
  def __init__ (self, chr, start, end, strand):
    self.chr   = chr
    self.start = int(start)
    self.end   = int(end)
    self.strand = strand
    self.meth = 0
  def printAll(self):
      attrs = vars(self)
      print ', '.join("%s: %s" % item for item in attrs.items())
  def matchMath(self, meth):
      if self.name == meth.name:
          self.meth = meth.val
      else:
          self.meth = "NaN"

class Meth:
    def __init__ (self, id, name, val, chr, pos):
        self.id = id
        self.name = name
        self.chr   = chr
        self.pos = int(pos)
        if isfloat(val) and float(val) >= 0 and float(val) <=1:
            self.val = val
        else:
            self.val = 'NaN'
    def printAll(self):
        attrs = vars(self)
        print ', '.join("%s: %s" % item for item in attrs.items())

###------------end function--------
##--------load all meth filenames 

outlogf = open(outlog,'w')
fnArray = []
with open(inpc) as inpf:
  for line in inpf.readlines():
    fnArray.append(line.strip())
    
nmethSamples = len(fnArray)
print "Meth sample number \t" + str(nmethSamples)

###----loading_Gene_Information----
geneList = []
with open(inpg) as inpgf:
  for line in inpgf.readlines():
      [identifer,chrom,start,end,strand] = line.strip().split("\t")
      geneList.append(identifer)
nGene = len(geneList)
outlogf.write("number_of_genes\t" + str(nGene) + "\n")
print "number of genes:\t" +str(nGene)

#####--------------------------------loading meth data. 
idfile = -1
outputH = open(out,'w')
for fn in fnArray:
  try:
    fntempf = open(fn)
    methValue = [0] * len(geneList)
    idfile = idfile + 1
    if idfile ==0:
        outputH.write("barcode"+"\t" + '\t'.join(geneList) + "\n")             
    for i, line in enumerate(fntempf):
      if i < nesp:
        pass
      else:
        [identifier,value,name,chrom, position] = line.strip().split("\t")
        if not ";" in name:
            if name in geneList:
               try: 
                   float(value)
                   methValue[geneList.index(name)] = value
               except ValueError:
                   pass 
            else:
                continue
        else:
            names = name.split(";")
            for idx in [geneList.index(g) for g in names if g in geneList] :
                methValue[idx] = value

    fntempf.close()  
    outputH.write(fn+"\t" + '\t'.join(map(str,methValue) )+ "\n")
  except IOError as e:
    print "I/O error({0}): {1}".format(e.errno, e.strerror)
    sys.exit()
outputH.close()
outlogf.write("number_of_meth_level3_files\t" + str(idfile + 1) + "\n")  
print "#----DONE"
