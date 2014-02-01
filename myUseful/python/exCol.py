#!/usr/bin/env python
#J.HE
#Desp:extract one column from a tab delimited file for large files, skip every headerline start with #
#Input: -n <2:column number> -f <filename> -o <output filename>
#Output: output file names
##TODO: stage 1 only handle one column extraction, stage2 handler multiple colum extraction


import os,getopt,sys

argv = sys.argv[1:]
usage = "Usage : \n python bakFile.py -n <> -f <inputfile> -o <output file >" 
try:
    opts,args = getopt.getopt(argv,"hr:b:s:t:")
except getopt.GetoptError:
    print usage
      sys.exit(2)

for opt, arg in opts:
  if opt == '-h':
    print usage
    sys.exit()
  elif opt in ("-n"):
    numCol     = arg
  elif opt in ("-f"):
    inputFile  = arg
  elif opt in ("-o"):
    outFile    = arg

print numCol + "\t" + inputFile + "\t" + outFile

outH = open(outFile, 'w')
with open(inputFile,15038) as inputH:
  for line in inputH.readlines():
    if 
