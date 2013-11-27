#!/usr/bin/python
# input: <file contains the input file names for tcga snp level2 data>
# output: <file including snps after filtering> <file: stats about the filtering>
#J.HE

import os
import math
import linecache
import sys, getopt

#------------test--------
#inp = "input_test"
#outp = "output_test"
#outpstat = outp + "_stat"
#------------test end ---
#get command line args
argv = sys.argv[1:]
inp = ''
outp = ''
try:
    opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
except getopt.GetoptError:
  print 'Usage:makeMat_SNP.py -i <input file:list of level2 snp file name.> -o <outputfile: output snp matrix file name> '
  sys.exit(2)
for opt, arg in opts:
  if opt == '-h':
    print 'Usage:makeMat_SNP.py -i <input file:list of level2 snp file name> -o <outputfile: output snp matrix file name>'
    sys.exit()
  elif opt in ("-i", "--ifile"):
    inp = arg
  elif opt in ("-o", "--ofile"):
    outp = arg
    outpstat = outp + "_stat"

conf_cut = 0.1
maf_cut = 0.05
fisher_cut =  10**(-8)

print "input file:" + inp
print "output file:" + outp
print "output stat file:" + outpstat
print "confidence cutoff:" + str(conf_cut)
print "MAF cutoff:" + str(maf_cut)
print "HWE cutoff:" + str(fisher_cut)

# get the number of lines in one file
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
def sysTime():
	time=os.system("date|awk '{print $4}'")
	print str(time)

def get_nth_line(f,n):
	# with open(f) as ff:
	line = linecache.getline(f,n)
	[key, gt, conf] = line.split("\t")
	if n % 10 == 0 and n > 10:
		linecache.clearcache()
	return key,gt,conf.strip()

nsample = file_len(inp)
sampleSize_cut = int(0.1 * nsample) 
with open(inp) as inf:
		for i, line in enumerate(inf):
			if i == 1:
				line_crt = line.strip()
				nsnp = file_len(line_crt)
				break

outp = outp +"_" + str(nsample) + ".mat"
outpstat = outpstat +"_" +  str(nsample) + ".txt"

def get_MAF(gtlist):
	return 1 - sum(map(int, gtlist)) / (2*len(gtlist))

def logfac(n):
	return math.log(math.factorial(int(n)))

def fisherTest(gtlist):
	N = len(gtlist)
	n_aa = len([i for i in range(len(gtlist)) if(gtlist[i] == 0)])
	n_ab = len([i for i in range(len(gtlist)) if(gtlist[i] == 1)])
	n_bb = len([i for i in range(len(gtlist)) if(gtlist[i] == 2)])
	prob = n_ab * math.log(2) + logfac(N) - \
		   logfac(n_aa) - logfac(n_ab) - logfac(n_bb) \
	 		+ logfac(2*n_aa + n_ab) + logfac(2 * n_bb + n_ab) - logfac(2 * N)
	return math.exp(prob)

cnt_f1 = 0
cnt_f2 = 0
cnt_f3 = 0
fout = open(outp,'w')
#the header
with open(inp) as inf:
  allsamples = [" "] * nsample
  for i, line in enumerate(inf):
  	line_crt = line.strip()
	allsamples[i] = line_crt
  fout.write("Identifier\t"+"\t".join(allsamples) + "\n")    


farray=[]
with open(inp) as f:
	for line in f.readlines():
		line_crt = line.strip()
		farray.append(line_crt)	

a=[open(i,'r',81920) for i in farray]
for k in range(nsnp):
	gt = [0] * nsample
	conf = [0] * nsample
	if k < 2:
		for f in a:
			line = f.next()
			pass
		continue
	else :
		ncol = 0
		for f in a:
			line = f.next()
			[key, g, c] = line.split("\t")
			gt[ncol] = int(g)
			conf[ncol] = float(c)
			ncol = ncol + 1
		samples=allsamples 
		allgt = gt		
		# 	filter1  confidence 
		temp = [j for j in range(len(conf)) if(conf[j] < conf_cut)]
		if len(temp) >= (nsample -sampleSize_cut):
		#  	filter2  maf
		   	maf = get_MAF(gt)
			if(maf > maf_cut):
		#	filter3 fisher test
				ftest = fisherTest(gt)
				if ftest > float(fisher_cut):
					fout.write(key +"\t"+ "\t".join(map(str,gt)) + "\n")
				else :
					cnt_f3 = cnt_f3 + 1
					continue
			else:
				cnt_f2 = cnt_f2 + 1
				continue
		else:
			cnt_f1 = cnt_f1 + 1
			continue



fout2 = open(outpstat,'w')
fout2.write("Number_of_total_snp" + "\t" + "filtered_out_Missing_GT" + \
	"\t" + "filtered_out_MAF" + \
	"\t"+ "filtered_out_HWE_Test" + \
	"\t" + "Number_of_final" + "\n")
fout2.write(str(nsnp) + "\t" + str(cnt_f1) + "\t" + str(cnt_f2) + "\t" + \
	str(cnt_f3) + "\t" + \
	str(nsnp - cnt_f1 - cnt_f2 - cnt_f3) + "\n")
	# str(nsnp - cnt_f1 - cnt_f2) + "\n")
fout.close()
fout2.close()


