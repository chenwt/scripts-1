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
import generalUtils as gu
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

def getCommElementIndex(list1,list2):
   ''' find the index of common elements for each list (keep the order \
           target list2 '''
   want = set(list2)
   idx1 =  [idx for (idx, val) in enumerate(list1) if val in want ]
   list1want = map(list1.__getitem__, idx1)
   idx2 =  [idx for (idx,val) in enumerate(want) if val in list1want ] 
   return [idx1, idx2]
   

def getIndiceforOneChrom(tarPlist, tempPlist):
    return(map(lambda x:getCommElementIndex(x[0], x[1]),\
               zip(tarPlist, tempPlist)) )
         

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

probelist = [[] for _ in range(25)]
genelist = [[] for _ in range(25)]
with open( tempglist ) as inpgf:
  for line in inpgf.readlines():
      probe, gene, chrom, ps  = line.strip().split("\t") 
      nChrom = gu.chr2Num(chrom) - 1
      probelist[nChrom].append(probe)
      genelist[nChrom].append(gene)

nprobe = reduce(lambda a,b: a+b, map(len,probelist))
ngene = reduce(lambda a,b: a+b, map(len,genelist))
outlogf.write("number_of_genes\t" + str(ngene) + "\n")
outlogf.write("number of uniq genes :\t" + str(nprobe) + "\n")
print "number of genes:\t" +str(ngene)
print "number of probes :\t" + str(nprobe) 

#####--------------------------------loading meth data. 
idfile = -1
outputH = open(out,'w')
idList = []
for fn in fnArray:
  try:
    f = open(fn)
    print fn
    methValue = [ [0] * i  for i  in map(len, probelist)] 
    tempIndices = [ [] for _ in range(25)] 
    tempProbelist  = [ [] for _ in range(25) ]
    tempVallist  = [ [] for _ in range(25) ]
    idfile = idfile + 1
    if idfile == 0:
        outputH.write("barcode"+"\t" + '\t'.join(\
                    reduce(list.__add__,probelist)) + "\n")             
        outputH.write("barcode"+"\t" + '\t'.join(\
                        reduce(list.__add__, genelist)) + "\n")             
    line = f.readline()
    while line:
        if re.match("^cg",line):
            [identifier, value, name,chrom, position] = line.strip().split("\t")
            nChrom = gu.chr2Num(chrom) - 1 
            tempProbelist[nChrom].append(identifier) 
            tempVallist[nChrom].append(value) 
        line = f.readline()
    f.close()  
    if re.findall("/", fn):
        sampleName = re.split("/", fn)[-1]
    else:
        sampleName = fn

    for chrIter in range(len(probelist)):
        tempProbelistIter = tempProbelist[chrIter] 
        tempVallistIter = tempVallist[chrIter]
        probelistIter = probelist[chrIter] 
        methValueIter = methValue[chrIter]
        for p in tempProbelistIter:
            if  p in probelistIter:
                methValueIter[probelistIter.index(p)] = \
                tempVallistIter[tempProbelistIter.index(p)] 
            else:
                continue 
    outputH.write( sampleName +"\t" +                                   \
                  '\t'.join(map(str,reduce(list.__add__,methValue) )) + "\n")

  except IOError as e:
    print "I/O error({0}): {1}".format(e.errno, e.strerror)
    sys.exit()
outputH.close()
outlogf.write("number_of_meth_level3_files\t" + str(idfile + 1) + "\n")  
print "#----DONE---"
