#!/usr/bin/env python
#J.HE
#Desp: used to get single start and end for each gene instead of transcripts for entrez_annotation_hg19.txt file
#input: entrez_annotation_hg19.txt
#output: entrez_annotation_hg19_unique.txt
##COMMENT: inputfile well sorted by chromosome

import sys,os
import getopt

inpFile = "entrez_annotation_hg19.txt"
outFile = "entrez_annotation_hg19_unique.txt"
argv = sys.argv[1:]
usage = 'python " + sys.argv[0] + " -i <input>  -o <output>'
example = 'python " + sys.argv[0] + " -i <input>  -o <output>'
try:
  opts,args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
except getopt.GetoptError:
  print usage + "\n" + example 
  sys.exit(2)
for opt, arg in opts:
  if opt == '-h':
    print usage + "\n" + example 
    sys.exit()
  elif opt in ("-i","--ifile"):
    inpFile = arg
  elif opt in ("-o","--ofile"):
    outFile = arg
    outlog = outFile + ".log"
print ('Script path'+ sys.argv[0])
print('Input file:' + inpFile)
print('Output file:'+ outFile)
print('Log file:'+ outlog)


outFileHandler = open(outFile,'w')
cntLine = 0 
with open(inpFile) as f:
  for i,line in enumerate(f):
    if i  == 0:
      geneInfoOld = ["gene","chrom","start","end","strand"]
    else:
      [geneId, chrom, start, end, strand] = line.strip().split("\t")
      # chrom = chrom.replace("chr","") 
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

