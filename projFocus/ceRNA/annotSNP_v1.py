#!/usr/bin/env python
#input:	<file: .mat file with snp/gene one each line to be annotated> 
#	<file: database file including the annotation infor: identifier,chr,pos,strand,dbsnpid>
#output:<file: annotat snp/gene with chr,start,end,strand,replace identifier with dbsnpid>
#TODO: ADD codes to do quality checking for database file! 
#examples: ~/tools/python/Python_current/python annot_SNP.py -i ~/SCRATCH/projFocus/ceRNA/data/snpArray/brca_snp_brca_snp_level3_839.mat -d ~/SCRATCH/database/projFocusRef/annot_GenomeWideSNP_6_5cols_clean.txt -o ~/SCRATCH/projFocus/ceRNA/data/snpArray/brca_snp_839.mat.annot  

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
outlogf = open(outlog,"w")

outlogf.write('Input file:' + inp1 + "\n")
outlogf.write('Database file:' + inp2 + "\n")
outlogf.write('Output file:'+ outp +"\n")
outlogf.write('Log file:'+ outlog + "\n")

outpf = open(outp,"w")
#load data file to be annotated
inp1Dict = {}
with open(inp1) as inf:
  for i,line in enumerate(inf):
    if i == 0:
      outHeader1 = line.split("\t",1)[1]
    elif i >0 :
      [key,val]= line.split("\t",1)
      inp1Dict[key] = val

outlogf.write("input snps:\t" + str(i) + "\n")
#load and annotate using database file
cnt_out = 0
outKey = []
with open(inp2) as inf2:
  for i, line in enumerate(inf2):
    if i == 0:
      outHeader2 = line.strip()
      # print outHeader2 + "\t" + outHeader1
      outpf.write(outHeader2 + "\t" + outHeader1)
    else :
      [key,info1,info2,info3,dbsnpid]  = line.strip().split("\t")
    # print key,info1,info2,info3
      if key in inp1Dict:
	cnt_out = cnt_out + 1
        outpf.write(dbsnpid + "\t" + info1 + "\t" + \
                    info2 + "\t" + info3+ "\t" + inp1Dict[key])
        # print i, key,info1,info2,info3
	outKey.append(key)
      else :
        continue
outlogf.write("drop snps:\t" + str(len(inp1Dict.keys()) - len(outKey)) + "\n")
outlogf.write("output snps:\t" + str(cnt_out) + "\n")

#for key in inp1Dict.keys() :
#  if key not in outKey:
#    outlogf.write("drop snp:\t" + key + "\n")

outpf.close()
outlogf.close()
