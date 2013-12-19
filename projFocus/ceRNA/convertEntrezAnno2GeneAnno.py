#!/usr/bin/env python
#J.HE.
#Desp: used to convert entrez id in entrez_annotation_hg19.txt to gene symbol with the same annotation
#input: <file: entrez_annotation_hg19.txt> <file: entrez2gene.txt>
#output: <file: gene_annotation_start_end_hg19.txt>

import os,sys
import getopt
usage = "python convertEntrezAnno2GeneAnnot.py -i <> -o <out.mat> -m <selected_sample_TCGA_barcode.txt> -o <?>"
argv = sys.argv[1:]
mapFile = ""
try:
      opts,args = getopt.getopt(argv,"hi:m:o:")
except getopt.GetoptError:
      print usage
      sys.exit(2)
for opt, arg in opts:
  if opt == "-h":
    print ""
  elif opt == "-i":
    inpFile = arg
  elif opt == "-m":
    mapFile = arg
  elif opt == "-o":
    outFile = arg

if mapFile == "":
  mapFile ="/ifs/scratch/c2b2/ac_lab/jh3283/database/geneIDconverter/entrez2gene.txt"

mapDict = {}

print ("loading mapping file...")
with open(mapFile) as f:
  for line in f.readlines():
    [oldId, newId] = line.strip().split("\t")
    mapDict[oldId] = newId

skipHeader = "TRUE"
outFileHandler = open(outFile,'w')
print("mapping ...")
cntUnmapped = 0
with open(inpFile) as f:
  for line in f.readlines():
    if skipHeader == "TRUE":
      skipHeader = "FALSE"
      outFileHandler.write(line)
    else:
      [oldId, value] = line.strip().split("\t",1)
      if mapDict.get(oldId,0):
	outFileHandler.write(mapDict.get(oldId) + "\t" + value + "\n")
      else:
	cntUnmapped = cntUnmapped + 1
	print oldId + "\t"

print " unmapped records number: \t" + str(cntUnmapped) 
print "#---END-----"  
outFileHandler.close()

