#!/usr/bin/env python
#input:<file: snp/gene one each line to be annotated> <file: database file including the annotation infor,chr,start,end,strand/if gene,use entrez id>
#output:<file: annotat snp/gene with chr,start,end,strand>
#TODO:

import os
import sys, getopt


argv = sys.argv[1:]
inp1 = ''
inp2 = ''
outp = ''
try:
  opts,args = getopt.getopt(argv,"hi:d:o:",["ifile=","dbfile'","ofile="])
except getopt.GetoptError:
  print 'python annot_SNP.py -i <inputfile> -d <snp annotation file> -o <output file name>' 
  sys.exit(2)

for opt, arg in opts:
  if opt == '-h':
    print 'python annot_SNP.py -i <inputfile> -d <snp annotation file> -o <output file name>' 
    sys.exit()
  elif opt in ("-i","--ifile"):
    inp1 = arg
  elif opt in ("-d","--dbfile"):
    inp2 = arg
  elif opt in ("-o","--ofile"):
    outp = arg
    outlog = outp + ".log"

print('Input file:' + inp1)
print('Database file:' + inp2)
print('Output file:'+ outp)
print('Log file:'+ outlog)
# ~/scripts/projFocus/ceRNA/test
# inp1 = '/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result_snp/run1/brca_snpGT_run1.mat'
# inp2 = '/ifs/scratch/c2b2/ac_lab/jh3283/database/affymetrix/GPL6801-4019.txt.bed'
# outp = inp1 + ".annot"
# outlog = outp + ".log"
outlogf = open(outlog,"w")

outlogf.write('Input file:' + inp1)
outlogf.write('Database file:' + inp2)
outlogf.write('Output file:'+ outp)
outlogf.write('Log file:'+ outlog)

outpf = open(outp,"w")

inp1Dict = {}
with open(inp1) as inf:
  for i,line in enumerate(inf):
    if i == 0:
      outHeader1 = line.split("\t",1)[1]
    elif i >0 :
      [key,val]= line.split("\t",1)
      inp1Dict[key] = val

with open(inp2) as inf2:
  for i, line in enumerate(inf2):
    if i == 0:
      outHeader2 = line.strip()
      # print outHeader2 + "\t" + outHeader1
      outpf.write(outHeader2 + "\t" + outHeader1)
    else :
      [key,info1,info2,info3]  = line.strip().split("\t")
    # print key,info1,info2,info3
      if key in inp1Dict:
        outpf.write(key + "\t" + info1 + "\t" + \
                    info2 + "\t" + info3+ "\t" + inp1Dict[key])
        # print i, key,info1,info2,info3
      else :
        continue

outpf.close()
outlogf.close()
