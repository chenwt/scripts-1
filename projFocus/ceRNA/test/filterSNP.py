#!/usr/bin/python
# input: <file contains the input file names>
# output: <file including snps after filtering>

import os
import math

# for one line
inp = "input_filter_GT.txt"
outp = "out_filter_GT.txt"
conf_cut = 0.01
maf_cut = 0.05
fisher_cut =  0.001
# sampleSize_cut = 30 # in real
sampleSize_cut = 2 # for test

# get the number of lines in one file
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

# line_crt = "test3.birdseed"
# i = 0
# get the nth line (snpid,GT,call_confidence) of one file
def get_nth_line(f,n):
	with open(f) as ff:
		n = int(n)
		for i, line in enumerate(ff):
			if i == n:
				[key, gt, conf] = line.split("\t")
	return key,gt,conf.strip()
	# ff.close()
	

# nrow = os.system("wc -l")
# nrow = int(nrow.strip())

nsample = file_len(inp)
with open(inp) as inf:
		for i, line in enumerate(inf):
			if i == 1:
				line_crt = line.strip()
				nsnp = file_len(line_crt)
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
for k in range(2, nsnp):
	if (k % 10000) == 0:
		print "1w SNPs processed..."
	gt = [0] * nsample
	conf = [0] * nsample
	samples = [" "] * nsample
	# extract information of one SNP across all samples
	with open(inp) as inf:
		for i, line in enumerate(inf):
			line_crt = line.strip()
			[key,g,c] = get_nth_line(line_crt,k)
			# print g + c 
			samples[i] = line_crt
			gt[i] = int(g)
			conf[i] = float(c)	
	if k == 2:
		# print "\t".join(samples)
		fout.write("\t".join(samples) + "\n")
	allsamples = samples
	allgt = gt		
# 	filter  confidence 
	a = [j for j in range(len(conf)) if(conf[j] < conf_cut)]
	if len(a) > sampleSize_cut:
		# print [ gt[i] for i in a ]
		# print [ samples[i] for i in a ]
		gt = [ gt[i] for i in a ]
		samples = [ samples[i] for i in a ]
#  	filter  maf
	   	maf = get_MAF(gt)
		if(maf > maf_cut):
#	filter fisher test
			ftest = fisherTest(gt)
			if ftest > fisher_cut:
				# print " ".join(map(str,gt))
				# print " ".join(samples)
				a = [j for j in range(len(allsamples)) if allsamples[j] not in samples]
				for i in a:
					allgt[i] = -2 
				# print key + "\t".join(map(str,allgt))
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
print "Number_of_snp" + "\t" + "filtered_out1" + "\t" + "filtered_out2" + "\t"+ "filtered_out3"
print nsnp,cnt_f1,cnt_f2,cnt_f3
fout.close()

	# print " ".join(map(str, gt))
	# print " ".join(map(str, conf))
	# k = k + 1
	# i = 0 

