#!/usr/bin/python
#J.HE
##input: <bam file full path> <region length>
##output: <split region file> 

import os,getopt
import sys

argv = sys.argv[1:]
input = ''
output = ''
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
    input = arg
  elif opt in ("-o","--ofile"):
    output = arg
    outlog = output + ".log"
print ('Script path'+ sys.argv[0])
print('Input file:' + input)
print('Output file:'+ output)
print('Log file:'+ outlog)

sys("samtools view -H bamfile")
