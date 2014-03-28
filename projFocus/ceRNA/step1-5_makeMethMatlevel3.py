#!/usr/bin/env python
# -*- coding: utf-8 -*-
#input: 1. <file:file names of the level 3 meth file> 2. <file: genes,and location of interest:format: identifier: chr. pos. strand>
#output: 1. <file: .mat.anno file including all meth info for the input files>
#Sample  Chromosome      Start   End     Num_Probes      Segment_Mean
#BITES_p_TCGAb61_SNP_S_GenomeWideSNP_6_C01_697154        1       61735   82650   11      0.5307
#COMMENT: debugged with barcode name;changed to include all probes(not targetd
# genes) 
      

import os, re
import sys, getopt
from subprocess import Popen, PIPE

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


##-----setting parameters
nesp	   = 2  ##number of lines to escape at the begining of each meth file
###-----functions------
def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False
    
###------------end function--------
##--------load all meth filenames 

outlogf = open(outlog,'w')
fnArray = []
with open(inpc) as inpf:
  for line in inpf.readlines():
    fnArray.append(line.strip())
    
nmethSamples = len(fnArray)
print "Meth sample number \t" + str(nmethSamples)
outlogf.write("Meth sample number \t" + str(nmethSamples) +"\n")

###----loading_probe_Information----
tempglist =  inpg + ".temp" 
cmd = "~/tools/python/Python_current/python \
/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/getGeneMethlevel3Probe.py \
-g " + inpg + " -p \
/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/03102014/methy/tumor/jhu-usc.edu_BRCA.HumanMethylation27.3.lvl-3.TCGA-BH-A0B0-01A-21D-A112-05.txt \
-o " + tempglist  
p = Popen(cmd, stdout = PIPE, stderr = PIPE, shell = True)
err = p.communicate()[1] 

if err:
    print "ERR in get cg list"
    sys.exit()
else:
    print "get probe success"

geneList = []
genelistreal = []
with open( tempglist ) as inpgf:
  for line in inpgf.readlines():
      probe, gene  = line.strip().split("\t")
      geneList.append(probe)
      genelistreal.append(gene)
nGene = len(geneList)
outlogf.write("number_of_genes\t" + str(nGene) + "\n")
outlogf.write("number of uniq genes :\t" + str(len(genelistreal)) + "\n")
print "number of probes:\t" +str(nGene)
print "number of uniq genes :\t" + str(len(set(genelistreal))) 

#####--------------------------------loading meth data. 
idfile = -1
outputH = open(out,'w')
idList = []
for fn in fnArray:
  try:
    fntempf = open(fn)
    print fn
    methValue = [0] * len(geneList)
    idfile = idfile + 1
    if idfile ==0:
        outputH.write("barcode"+"\t" + '\t'.join(geneList) + "\n")             
        outputH.write("barcode"+"\t" + '\t'.join(genelistreal) + "\n")             
    for i, line in enumerate(fntempf):
      if i < nesp:
        pass
      else:
        [identifier, value, name,chrom, position] = line.strip().split("\t")
        if identifier in geneList:
               try: 
                   value=float(value)
                   temp = methValue[geneList.index(identifier)] 
                   methValue[geneList.index(identifier)] = value         
               except ValueError:
                   pass 
        else:
            continue
    fntempf.close()  
    if re.findall("/", fn):
        sampleName = re.split("/", fn)[-1]
    else:
        sampleName = fn
    outputH.write( sampleName +"\t" + '\t'.join(map(str,methValue) )+ "\n")
  except IOError as e:
    print "I/O error({0}): {1}".format(e.errno, e.strerror)
    sys.exit()
outputH.close()
outlogf.write("number_of_meth_level3_files\t" + str(idfile + 1) + "\n")  
print "#----DONE"
