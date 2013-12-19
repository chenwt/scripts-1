#!/usr/bin/env python
#J.HE
#Desp: used to get single start and end for each gene instead of transcripts for entrez_annotation_hg19.txt file
#input: entrez_annotation_hg19.txt
#output: entrez_annotation_hg19_unique.txt
##COMMENT: inputfile well sorted by chromosome

import sys
inpFile = "entrez_annotation_hg19.txt"
outFile = "entrez_annotation_hg19_unique.txt"

outFileHandler = open(outFile,'w')
cntLine = 0 
with open(inpFile) as f:
  for i,line in enumerate(f):
    if i  == 0:
      geneInfoOld = ["gene","chrom","start","end","strand"]
    else:
      [geneId, chrom, start, end, val2, strand, val3] = line.strip().split("\t",6)
      chrom = chrom.replace("chr","") 
      geneInfoCrt = [geneId, chrom, start, end, strand]
      if not geneInfoOld: 
        geneInfoOld = geneInfoCrt 
	continue
      elif not geneInfoOld and geneInfoOld[0] == geneId : 
        if geneInfoCrt[4] == "-":
          if int(geneInfoCrt[2]) > int(geneInfoOld[2]):
            geneInfoOld[2]  = start
	  elif int(geneInfoCrt[3]) < int(geneInfoOld[3]) :
            geneInfoOld[3]  = end
	  continue
        else :
          if int(geneInfoCrt[2]) < int(geneInfoOld[2]) :
            geneInfoOld[2]  = start
          elif int(geneInfoCrt[3]) > int(geneInfoOld[3]) :
            geneInfoOld[3]  = end
	  continue
      elif geneInfoOld and not geneInfoOld[0] == geneId : 
        # compare to update or not 
        outFileHandler.write( "\t".join(geneInfoOld) + "\n" )
	geneInfoOld = geneInfoCrt 
    cntLine = cntLine + 1

print "line processed:\t" + str(cntLine)
print "#----END---"
outFileHandler.close()

