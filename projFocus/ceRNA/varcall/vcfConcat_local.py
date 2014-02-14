#!/usr/bin/python
#J.HE

import sys,getopt
import os
import re

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
print ('Script path'+ sys.argv[0])
print('Input file:' + input)
print('Output file:'+ output)

##input is a pattern of files in the current dir
dir ='.'
pattern = input
files = [f for f in os.listdir(dir) if os.path.isfile(os.path.join(dir,f)) and re.findall(pattern,f) ] 
cntf  = 0

outputH = open(output,'w')
for f in files:
  print "file:\t" + f
  with open(f) as fH:
    for line in fH.readlines():
      if cntf == 0 and re.findall(r'^#',line):
	outputH.write(line)
      else :
	if re.findall(r'^#',line):
	  pass
	else:
	  outputH.write(line)
  cntf = cntf + 1
 
outputH.close()
print "##------END---"


