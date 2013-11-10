#!/usr/bin/python
#input: <meth27.mat> <meth460.mat>
#ouput: <meth.mat>

import os
import sys,getopt

argv = sys.argv[1:]
inp1 = ''
inp2 = ''
try:
  opts,args = getopt.getopt(argv,"hi:m:o:",["meth27=","meth450='","out="])
except getopt.GetoptError:
  print 'python compareFileCols.py -i <meth27 file> -m <meth450 file> -o <out file name>'
  sys.exit(2)
for opt, arg in opts:
  if opt == '-h':
    print 'python compareFileCols.py -i <meth27 file> -m <meth450 file> -o <out file name>'
    sys.exit()
  elif opt in ("-i","--meth27"):
    inp1 = arg
  elif opt in ("-o","--out"):
    outp = arg
    outlog = arg + ".log"
  elif opt in ("-m","--meth450"): 
    inp2 = arg

if inp1 == '' or inp2 == '' or outp == '':
  print "input parameter ERROR!"
  print "python compareFileCols.py -i <meth27 file> -m <meth450 file> -o <out file name> "
  sys.exit(2)
###----test_start---------
#inp1 = 'input_test_meth27.txt'
#inp2 = 'input_test_meth450.txt'
#outp = 'output_test_meth.mat'
#outlog = 'output_test_meth.mat.log'
###----test_end------------
outlogf = open(outlog, 'w')
inp1Dict = {}
with open(inp1) as inpf1:
  for i,line in enumerate(inpf1):
    [key,val] = line.strip().split("\t",1)
    inp1Dict[key] = val
  outlogf.write("Lines in " + inp1 + ":" + str(i + 1) + "\n" )

outpf = open(outp,'w')
count = 0 
with open(inp2) as inpf2:
  for i, line in enumerate(inpf2):
    [key,val] = line.strip().split("\t",1)
    if key in inp1Dict.keys():
      outpf.write(key + "\t" + inp1Dict[key] + "\t" + val + "\n")
      count = count + 1
    else:
      continue
  outlogf.write("Lines in " + inp2 + ":" + str(i + 1) + "\n" )

outlogf.write("Lines in " + outp + ":" + str(count) + "\n" )
outpf.close()
outlogf.close()
print "#-----DONE------"


