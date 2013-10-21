#!/usr/bin/python
# input: <file contains the input file names>
# output: <file including snps after filtering> <file: stats about the filtering>
#J.HE

import os
import math
import linecache

# for one line
# inp = "input_test.txt"
# inp = "input_filter_GT_top400.txt"
inp = "input_filter_GT_bottom436.txt"
outp = "out_filter_GT"
outpstat = "out_filter_GT_stat"
conf_cut = 0.01
maf_cut = 0.05
fisher_cut =  0.001
sampleSize_cut = 30 # in real
# sampleSize_cut = 2 # for test

	

# get the number of lines in one file
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
def sysTime():
	time=os.system("date|awk '{print $4}'")
	print str(time)
# line_crt = "test3.birdseed"
# i = 0
# get the nth line (snpid,GT,call_confidence) of one file
# def get_nth_line(f,n):
# 	with open(f) as ff:
# 		n = int(n)
# 		for i, line in enumerate(ff):
# 			if i == n:
# 				[key, gt, conf] = line.split("\t")
# 	return key,gt,conf.strip()

def get_nth_line(f,n):
	# with open(f) as ff:
	line = linecache.getline(f,n)
	[key, gt, conf] = line.split("\t")
	if n % 10 == 0 and n > 10:
		linecache.clearcache()
	return key,gt,conf.strip()

nsample = file_len(inp)
with open(inp) as inf:
		for i, line in enumerate(inf):
			if i == 1:
				line_crt = line.strip()
				nsnp = file_len(line_crt)
				break

outp = outp + str(nsample) + ".txt"
outpstat = outpstat + str(nsample) + ".txt"

def get_MAF(gtlist):
	return 1 - sum(map(int, gtlist)) / (2*len(gtlist))

def logfac(n):
	return math.log(math.factorial(int(n)))

def fisherTest(gtlist):
	N = len(gtlist)
	n_aa = len([i for i in range(len(gtlist)) if(gtlist[i] == 0)])
	n_ab = len([i for i in range(len(gtlist)) if(gtlist[i] == 1)])
	n_bb = len([i for i in range(len(gtlist)) if(gtlist[i] == 2)])
	# prob = (2 ** n_ab) * math.factorial(N) / \
	 # math.factorial(n_aa) / math.factorial(n_ab) / \
	  # math.factorial(n_bb) * math.factorial(2*n_aa + n_ab) * \
	   # math.factorial(2 * n_bb + n_ab) / math.factorial(2 * N)
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
  fout.write("\t".join(allsamples) + "\n")    

# sysTime()

farray=[]
with open(inp) as f:
	for line in f.readlines():
		line_crt = line.strip()
		farray.append(line_crt)	
# print " ".join(farray)
# for k in range(3, nsnp):
# 	if k % 100 == 0  and k > 100:
# 		print str(k) + " SNPs processed..."
# 	# sysTime()
# 	# extract information of one SNP across all samples
# 	with open(inp) as inf:
# 		for i, line in enumerate(inf):
# 			line_crt = line.strip()
# 			[key,g,c] = get_nth_line(line_crt,k)
a=[open(i,'r',81920) for i in farray]
for k in range(nsnp):
	gt = [0] * nsample
	conf = [0] * nsample
	if k < 2:
		for f in a:
			line = f.next()
			# print line
			pass
		continue
	else :
		ncol = 0
		for f in a:
			line = f.next()
			#samples[i] = line_crt
			[key, g, c] = line.split("\t")
			gt[ncol] = int(g)
			conf[ncol] = float(c)
			ncol = ncol + 1
		samples=allsamples 
		allgt = gt		
		# 	filter1  confidence 
		temp = [j for j in range(len(conf)) if(conf[j] < conf_cut)]
		if len(temp) > sampleSize_cut:
			gt = [ gt[i] for i in temp ]
			samples = [ samples[i] for i in temp ]
		#  	filter2  maf
		   	maf = get_MAF(gt)
			if(maf > maf_cut):
		#	filter3 fisher test
				ftest = fisherTest(gt)
				if ftest > fisher_cut:
					temp2 = [j for j in range(len(allsamples)) if allsamples[j] not in samples]
					for i in temp2:
						allgt[i] = -2 
					fout.write(key + "\t".join(map(str,allgt)) + "\n")
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
fout2.write("Number_of_total_snp" + "\t" + "filtered_out1" + \
	"\t" + "filtered_out2" + \
	"\t"+ "filtered_out3" + \
	"\t" + "Number_of_final" + "\n")
fout2.write(str(nsnp) + "\t" + str(cnt_f1) + "\t" + str(cnt_f2) + "\t" + \
	str(cnt_f3) + "\t" + \
	str(nsnp - cnt_f1 - cnt_f2 - cnt_f3) + "\n")
	# str(nsnp - cnt_f1 - cnt_f2) + "\n")
fout.close()
fout2.close()


