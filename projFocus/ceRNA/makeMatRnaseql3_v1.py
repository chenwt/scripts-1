#!/usr/bin/python
# input:   <file: file name per line> 
#	  <string: the col number for data:default first column is the key,and every file has the same row number>
# output: <.mat file integrate> <file: stats about the filtering>
#J.HE

import os
import linecache
import sys, getopt

######## test use
## for one line
#inp = "input_test.txt"
#outp = "output_test.txt"
#outlog = "output_test.log"
######### test use end

#get command line args
argv = sys.argv[1:]
inp = ''
outp = ''
try:
  opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
except getopt.GetoptError:
  print 'test_rnaExpMat.py -i <inputfile> -o <outputfile>'
  sys.exit(2)
for opt, arg in opts:
  if opt == '-h':
     print 'test_rnaExpMat.py -i <inputfile> -o <outputfile>'
     sys.exit()
  elif opt in ("-i", "--ifile"):
     inp = arg
  elif opt in ("-o", "--ofile"):
     outp = arg
     outlog = outp + ".log"
print 'Input file is :', (inp)
print 'Output file is :', outp
print 'Log file is: ', outlog


def sysTime(fdlog):
	time=os.system("date|awk '{print $4}'")
	print str(time)
	fdlog.write(str(time))

# get the number of lines in one file
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

#get dimension of output matrix
cntCol = file_len(inp)
with open(inp) as inf:
		for i, line in enumerate(inf):
			if i == 1:
				line_crt = line.strip()
				cntRow = file_len(line_crt)
				break

fout = open(outp,'w')
#write  header
with open(inp) as inf:
  allsamples = [" "] * (cntCol + 1)
  allsamples[0] = "GeneSymbol"
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
# AATF|26574      4273.00 5.57160609223946e-05    uc002hni.2,uc002hnj.2
a=[open(i,'r',81920) for i in farray]
for k in range(cntRow):
	expt=[0] * cntCol
	if k < 1:
		for f in a:
			line = f.next()
			pass
		continue
	else :
		ncol = 0
		for f in a:
			line = f.next()
			[key, e] = line.split("\t")[:2]
			expt[ncol] = int(e)
			ncol = ncol + 1
			keynew = key.split("|")[0]
		if not keynew == '?':
		  fout.write(keynew + "\t"+"\t".join(map(str,expt)) + "\n")
		else:
		  continue

print "Number_of_gene" + "\t" + "Number_of_samples" + "\n"
print str(cntRow) + "\t" + str(cntCol) + "\n"
fout2 = open(outlog,'w')
fout2.write("Number_of_gene" + "\t" + "Number_of_samples" + "\n")
fout2.write(str(cntRow) + "\t" + str(cntCol) + "\n")
fout.close()
fout2.close()

print "-------DONE---------------"



