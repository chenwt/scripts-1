#!/usr/bin/env python
#input: <file: exp.mat>
#ouput: <file: exp.mat with single tss for each gene>

import sys, getopt
import os
argv = sys.argv[1:]

inpe = argv[0]
out = argv[1]
outf = open(out,'w')
geneDict = {}
cnt = 0
with open(inpe) as inpef:
  for i, line in enumerate(inpef):
    if i == 0:
      outf.write(line.strip())
    else:
      [identifier, chrom, pos, strand, val] = line.strip().split("\t", 4)
      if identifier in geneDict.keys():
	if (strand == "-" and int(pos) > int(geneDict[identifier][1]) ) or (strand == "+" and int(pos) < int(geneDict[identifier][1]) ):
	  cnt = cnt + 1
	  geneDict[identifier] = [chrom,pos,strand,val] 
	else:
	  pass
      else:
	geneDict[identifier] = [chrom,pos,strand,val]

for key in geneDict.keys():
  outf.write( key + "\t"  + "\t".join(geneDict[key]) + "\n")
print "number of non-First-exon deleted:\t" + str(cnt)
outf.close()
print "#------------END-----"
