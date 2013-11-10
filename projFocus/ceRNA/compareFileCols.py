#!/usr/bin/python
#input: <snp.anno> <exp.anno>
#ouput: <snp.ann with comm cols> <exp.anno with comm cols>

import os
import sys,getopt

argv = sys.argv[1:]
inp1 = ''
inp2 = ''
try:
    opts,args = getopt.getopt(argv,"hs:e:",["snpfile=","expfile'"])
except getopt.GetoptError:
    print 'python compareFileCols.py -s <snp annotation file> -e <exp annotation file> '
    sys.exit(2)
for opt, arg in opts:
  if opt == '-h':
    print 'python compareFileCols.py -s <snp annotation file> -e <exp annotation file> '
    sys.exit()
  elif opt in ("-s","--snpfile"):
    inp1 = arg
    outp1 = inp1 + ".comm"
  elif opt in ("-e","--expfile"): 
    inp2 = arg
    outp2 = inp2 + ".comm"

if inp1 == '' or inp2 == '':
  print "input parameter ERROR!"
  print "python compareFileCols.py -s <snp annotation file> -e <exp annotation file> "
  sys.exit(2)
###----test_start---------
#inp1 = 'input_test_snp.txt'
#inp2 = 'input_test_exp.txt'
#outp1 = 'output_test_snp.txt'
#outp2 = 'output_test_exp.txt'
###----test_end------------


with open(inp1) as inpf1:
  for i, line in enumerate(inpf1):
    if i == 0:
      colNames1 = line.strip().split("\t")[4:]
    else:
      break
with open(inp2) as inpf2:
  for i, line in enumerate(inpf2):
    if i == 0:
      colNames2 = line.strip().split("\t")[4:]
    else:
      break
#
def mySubstr(x):
  return(x[0:19])
colN1 = map(mySubstr,colNames1)
colN2 = map(mySubstr,colNames2)
comm = list(set(colN1)&set(colN2))
#print ",".join(comm)
idx_colN1 = [0,1,2,3]
idx_colN2 = [0,1,2,3]
for i in range(len(comm)):
  idx_colN1.append(colN1.index(comm[i]) + 4) 
for i in range(len(comm)):
  idx_colN2.append(colN2.index(comm[i]) + 4) 
#print ",".join(map(str,idx_colN2))

def outMyIdx(valArray,idxArray):
  if len(idxArray) < len(valArray):
    valOutArray = []
    for idx in idxArray:
      valOutArray.append(valArray[idx])
  return(valOutArray)
#-------test_start
#testval = [1,1,1,2,2,2,3,3,3]
#testidx = [0,1,2]
#print ":".join(map(str,outMyIdx(testval,testidx)))
#-------test_end

def outMyIdxFile(filename,out,idx):
  if len(idx_colN1) == len(idx_colN2):
    fout = open(out,'w')
    with open(filename) as f:
      for line in f.readlines():
	colsArray = line.strip().split("\t")
	colsArray_new = outMyIdx(colsArray,idx)
	#print "\t".join(colsArray_new)
	fout.write("\t".join(colsArray_new) + "\n")
  fout.close()

outMyIdxFile(inp1,outp1,idx_colN1)
outMyIdxFile(inp2,outp2,idx_colN2)

