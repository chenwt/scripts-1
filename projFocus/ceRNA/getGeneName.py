#!/usr/bin/env python
#J.HE
#input: <txt.file>
#output: <txt.file with uniqe gene names)

import os
import sys, getopt
import re

argv = sys.argv[1:]
inp1 = ''
inp2 = ''
outp = ''
try:
  opts,args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
except getopt.GetoptError:
  print 'python .py -i <inputfile>  -o <output file name>' 
  sys.exit(2)

for opt, arg in opts:
  if opt == '-h':
    print 'python annot_SNP.py -i <inputfile> -o <output file name>' 
    sys.exit()
  elif opt in ("-i","--ifile"):
    input = arg
  elif opt in ("-o","--ofile"):
    output = arg

outputH = open(output, 'w')
controllist = ['TCGA','DNA','RNA','SNP','LOH','COSMIC','II','GBM','COAD','IV','III','S1','S2','S3','S4','GISTIC','GST','GFP']
with open(input) as f:
  for line in f.readlines():
    templist = re.split('-|\s+|\t|\s|, |,|\. |\.|;|; |: |:|\'| \(|\) |/|\+|\-|\*',line.strip())
    res = [word for word in templist if (not re.match(r'\)\b',word)) and  (not re.match(r'^[0-9\(]',word)) and (not word in controllist) and  (not re.findall(r'\b\d+\b',word)) and len(word) > 1 and len(word) < 10 and re.findall(r"\b[A-Z0-9]+(?:\W*[A-Z0-9]+)*\b",word)]
    if len(res) > 0  :
      outputH.write("\n".join(list(set(res))) + "\n")
    else:
      continue 
    # outputH.write("\n".join([word for word in templist if len(word) > 1 and re.findall(r"\b[A-Z]+(?:\W*[A-Z]+)*\b",word)]) + "\n")
   
# sysCmd = "sort "+output+" | uniq >"+ output+".temp"
# os.system(sysCmd)
# # sysCmd = "rm "+ output
# # os.system(sysCmd)
# # sysCmd = "mv "+ output+".temp " + output
# # os.system(sysCmd)
print "#-----END---"
