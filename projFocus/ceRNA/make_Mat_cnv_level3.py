#!/usr/bin/env python
#input: <file:file names of the level 3 cnv file>
#output: 1. <file: .mat.anno file including all cnv info for the input files>
#Sample  Chromosome      Start   End     Num_Probes      Segment_Mean
#BITES_p_TCGAb61_SNP_S_GenomeWideSNP_6_C01_697154        1       61735   82650   11      0.5307
      

import os
import sys, getopt

argv = sys.argv[1:]
try:
  opts,args = getopt.getopt(argv, "h:i:o:")
except:
  print "ERROR:usage:-----"
  sys.exit()
for opt,arg in opts:
  if opt == '-h':
    print "ERROR:usage:-----"
    sys.exit()
  elif opt == '-i':
    inp = arg
  elif opt == '-o':
    out = arg
    outlog = out + ".log"

#------test_start----
#inp = "input_test_cnv_tu.txt"
#out = "output_test.txt"
#outlog = out + ".log"
#------test_end------
outlogf = open(outlog,'w')
nesp = 1
nval = 5 - 1
fnArray = []
with open(inp) as inpf:
  for line in inpf.readlines():
    fnArray.append(line.strip())

resDict = {}
idfile = -1
for fn in fnArray:
 try:
   fntempf = open(fn)
   idfile = idfile + 1
   for i, line in enumerate(fntempf):
     if i < nesp:
       pass
     else:
       [id,chr,pstart,pend,probenum,val] = line.strip().split("\t")
       key = str(chr) + "_" + str(pstart) + "_" + str(pend) 
       val = val + "_" + str(idfile)
       if key in resDict.keys():
	 resDict[key] = resDict[key] + [val]
       else :
	 resDict[key] = [val]
   fntempf.close()
   outlogf.write(str(i)+ fn + " processed...\n")
   print str(i) + fn + " processed..."
 except IOError as e:
       print "I/O error({0}): {1}".format(e.errno, e.strerror)
       sys.exit()

def mySplit_1(x):
  return(str(x).split("_")[0])
def mySplit_2(x):
  return(int(str(x).split("_")[1]))

outf = open(out,'w')
cnt_out = 0 
for key in resDict.keys():
  tempRow = [0] * (idfile + 1)
  for val in resDict[key]:
    tempRow[mySplit_2(val)]=mySplit_1(val)      
  outf.write(key + "\t" + "\t".join(map(str,tempRow)) + "\n")
  cnt_out = cnt_out + 1 

outlogf.write("Files processed: " + str(len(fnArray))) 
outlogf.write("Total cnv segs : " + str(cnt_out) )
print "#--------END------"



