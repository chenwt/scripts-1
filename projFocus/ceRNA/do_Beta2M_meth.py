#!/usr/bin/env python
#input: <file: meth.mat mergeing 27k and 450k, and annotation>
#output: <file: transform the original beta value to M value, and delete NA value(-10 as denoted in the original beta value>
#J.HE

import sys,getopt
import os
import math

inp1 = sys.argv[1]
####------test_start-------------
#inp1 = 'input_test_meth_all.txt'

#--------test_end-------------
outp = inp1 + ".Mval"
outlog = inp1 + ".Mval.log"

print "input file: " + inp1
print "output file: " + outp

def beta2M(x):
  x = float(x)
  return math.log(x / (1 - x), 2)

def flag(x):
  return(int(float(x)>0))

outpf = open(outp,'w')
outlogf = open(outlog,'w')
cnt_out = 1
with open(inp1) as inp1f:
  for i,line in enumerate(inp1f):
    if i == 0:
      outpf.write(line)
    else:
      line_temp = []
      line_temp = line.strip().split("\t")
      vals_b = line_temp[4:] 
      flags = map(flag,vals_b)
      if 0 not in flags:
	vals_m = map(beta2M,vals_b) 
	outpf.write("\t".join(line_temp[:4]) + "\t" + "\t".join(map(str,vals_m)) + "\n")
	cnt_out = cnt_out + 1
      else:
	continue
outlogf.write("Totoal input lines: " + str(i + 1) + "\n")
outlogf.write("Totoal output lines: " + str(cnt_out + 1) + "\n")
