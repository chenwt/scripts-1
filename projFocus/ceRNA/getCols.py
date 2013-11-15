#!/usr/bin/python
#input: <file: matfile> <file:colnames to retain>
#ouput: <file: matfile with specified column names> 
#COMMENT: the first 4 columnas are alwasy retained

import os
import sys,getopt

argv = sys.argv[1:]
inp1 = ''
inpc = ''
out = ''
try:
  opts,args = getopt.getopt(argv,"hi:c:o:")
except getopt.GetoptError:
  print 'python getCols.py -i <snp annotation file> -c <file:give the columns> -o <exp annotation file> '
  sys.exit(2)
for opt, arg in opts:
  if opt == '-h':
    print 'python getCols.py -i <snp annotation file> -c <file:give the columns> -o <exp annotation file> '
    sys.exit()
  elif opt in ("-i"):
    inp1 = arg
  elif opt in ("-c"):
    inpc = arg 
  elif opt in ("-o"): 
    out = arg

if inp1 == '' or inpc == '':
  print "input parameter ERROR!"
  print 'python getCols.py -i <snp annotation file> -c <file:give the columns> -o <exp annotation file> '
  sys.exit(2)
#------------loading file
colnames_new = []
with open(inpc) as inpcf:
  for line in inpcf.readlines():
      colnames_new.append(line.strip())
outf = open(out,'w')  
with open(inp1) as inp1f:
  for i, line in enumerate(inp1f):
    if i == 0:
      colnames_old = line.strip().split("\t")[4:]
      #index_col2keep = map(lambda x:x+4,[a for a in range(len(colnames_old)) if colnames_old[a] in colnames_new])
      index_col2keep = [0,1,2,3] + map(lambda x:x+4,[a for a in range(len(colnames_old)) if colnames_old[a] in colnames_new])
      #print "columns_to_keep\t" +  "\t".join(map(str,index_col2keep)) + "\n"
      colnames_old = line.strip().split("\t")
      #index_col2keep = map(lambda x:x+4,[a for a in range(len(colnames_old)) if colnames_old[a] in colnames_new])
      outf.write("\t".join([colnames_old[i] for i in index_col2keep]) + "\n")
    else:
      col_old = line.strip().split("\t")
      col_new = [col_old[i] for i in index_col2keep]  
      outf.write("\t".join(map(str,col_new)) + "\n")
 

outf.close()
print str(len(colnames_new)) + " columns kept..."
print "#-----END-----"

