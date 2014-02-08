#!/usr/bin/python
# input:   <file: file name per line> 
#	  <string: the col number for data:default first column is the key,and every file has the same row number>
# output: <.mat file integrate> <file: stats about the filtering>
#TODO: MAKE SURE THE MEMORY SIZE IS ENOUGTH BEFORE RUNNING THIS
#J.HE

import os
import re 
import linecache
import sys, getopt

######## test use
## for one line
#inp = "input_test.txt"
#outp = "output_test.txt"
#outlog = "output_test.log"
#nval = 2
######### test use end

#get command line args
argv = sys.argv[1:]
inp = ''
outp = ''
nval = 2
usage = 'python sys.argv[1] ' + ' -i <inputfile> -g <annotated gene file> -c <2:column of value> -o <outputfile> -e <2: rows to escape for each file>'
try:
  opts, args = getopt.getopt(argv,"hi:c:e:o:",["ifile=","valCol=","rowEsp=","ofile="])
except getopt.GetoptError:
  print usage
  sys.exit(2)
for opt, arg in opts:
  if opt == '-h':
    print
    sys.exit()
  elif opt in ("-i", "--ifile"):
     inp = arg
  elif opt in ("-c", "--ifile"):
     nval = int(arg)
  elif opt in ("-g", "--ifile"):
     inputgene = int(arg)
  elif opt in ("-e", "--ifile"):
     nesp = int(arg)
  elif opt in ("-o", "--ofile"):
     outp = arg
     outlog = outp + ".log"
print 'Input file is :', inp
print 'Output file is :', outp
print 'Log file is: ', outlog


# get the number of lines in one file
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

#get dimension of output matrix
cntCol = file_len(inp)
with open(inputgene) as f:
  for line in f.readlins()
  # assume no header for input gene file
    genelist.append(line.strip().split("\t",1)[0]
cntRow = len(genelist)


fout = open(outp,'w')
#write  header
with open(inp) as inf:
  allsamples = [" "] * (cntCol + 1)
  allsamples[0] = "gene"
  for i, line in enumerate(inf):
  	line_crt = line.strip()
	allsamples[(i+1)] = line_crt
  fout.write("\t".join(allsamples) + "\n")    

#filehandler array 
farray=[]
with open(inp) as f:
	for line in f.readlines():
		line_crt = line.strip()
		farray.append(line_crt)	

# data
#81920 for 400 file,163840 
a=[open(i,'r',163840) for i in farray]
for k in range(cntRow):
	expt=[0] * cntCol
	if k < nesp:
		for f in a:
			line = f.next()
			pass
		continue
	else :
		ncol = 0
		for f in a:
			line = f.next()
			keynew = line.split("\t")[2] 
			if not ";" in keynew and keynew in genelist: 
			  
			val = line.split("\t")[(nval - 1)]
			regexp = r'[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?'
			if re.match(regexp,val):
			  expt[ncol] = float(val)
			else :
			  expt[ncol] = "NaN"
			  #remember this is to denote missing value!
			ncol = ncol + 1
		fout.write(keynew + "\t"+"\t".join(map(str,expt)) + "\n")

print "Number_of_Identifier" + "\t" + "Number_of_samples" + "\n"
print str(cntRow - nesp + 1) + "\t" + str(cntCol) + "\n"
fout2 = open(outlog,'w')
fout2.write("Number_of_Identifier" + "\t" + "Number_of_samples" + "\n")
fout2.write(str(cntRow - nesp + 1) + "\t" + str(cntCol) + "\n")
fout.close()
fout2.close()

print "-------DONE---------------"
