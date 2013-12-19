#!/usr/bin/env python
#input:	<file: .mat file with snp/gene one each line to be annotated> 
#	<file: database file including the annotation infor: identifier,chr,pos,end,strand>
#output:<file: annotat snp/gene with chr,start,end,strand>
#TODO: ADD codes to do quality checking for database file! 
#examples: ~/tools/python/Python_current/python annotSomMutation.py -i ~/SCRATCH/projFocus/ceRNA/data/somatic/somtic.mat -d ~/SCRATCH/database/projFocusRef/gene_annotation_hg19_unique_start_end.txt -o ~/SCRATCH/projFocus/ceRNA/data/snpArray/brca_snp_839.mat.annot  

import os
import sys, getopt


argv = sys.argv[1:]
inp1 = ''
inp2 = ''
outp = ''
usage = 'python annot_SNP.py -i <inputfile> -d <snp annotation file> -o <output file name>' 

try:
  opts,args = getopt.getopt(argv,"hi:d:o:",["ifile=","dbfile'","ofile="])
except getopt.GetoptError:
  print usage 
  sys.exit(2)

for opt, arg in opts:
  if opt == '-h':
    print usage 
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
      print outHeader1
    else :
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
      print outHeader2 + "\t" + outHeader1
      outpf.write(outHeader2 + "\t" + outHeader1)
    else :
      [key,info1,info2,info3,info4]  = line.strip().split("\t")
    # print key,info1,info2,info3
      if key in inp1Dict:
	cnt_out = cnt_out + 1
        outpf.write(key + "\t" + info1 + "\t" + \
                    info2 + "\t" + info3+ "\t" + info4 + "\t" + inp1Dict[key])
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
