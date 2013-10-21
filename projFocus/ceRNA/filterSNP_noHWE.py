#!/usr/bin/python
# input: <file contains the input file names>
# output: <file including snps after filtering> <file: stats about the filtering>

import os
import math

# for one line
inp = "input_filter_GT_top400.txt"
#inp = "input_test"
outp = "out_filter_GT"
outpstat = "out_filter_GT_stat"
conf_cut = 0.01
maf_cut = 0.05
#fisher_cut =  0
sampleSize_cut = 30 # in real
# sampleSize_cut = 2 # for test

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

nsample = file_len(inp)
with open(inp) as inf:
		for i, line in enumerate(inf):
			if i == 1:
				line_crt = line.strip()
				nsnp = file_len(line_crt)
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
time1=os.system("date|awk '{print $4}'")
print "time1" + str(time1)
with open(inp) as inf:
  allsamples = [" "] * nsample
  for i, line in enumerate(inf):
	line_crt = line.strip()
	allsamples[i] = line_crt
  fout.write("\t".join(allsamples) + "\n")    

time2=os.system("date|awk '{print $4}'")
print "time2" + str(time2)
for k in range(2, nsnp):
	#if k % 100 == 0  and k > 100:
	#	print str(k) + " SNPs processed..."
	gt = [0] * nsample
	conf = [0] * nsample
	time3=os.system("date|awk '{print $4}'")
	print "time3" + str(time3)
	#samples = [" "] * nsample
	# extract information of one SNP across all samples
	with open(inp) as inf:
		for i, line in enumerate(inf):
			line_crt = line.strip()
			[key,g,c] = get_nth_line(line_crt,k)
			# print g + c 
			#samples[i] = line_crt
			gt[i] = int(g)
			conf[i] = float(c)	
	#if k == 2:
	#	fout.write("\t".join(samples) + "\n")
	samples=allsamples 
	allgt = gt		
# 	filter  confidence 
	a = [j for j in range(len(conf)) if(conf[j] < conf_cut)]
	if len(a) > sampleSize_cut:
		gt = [ gt[i] for i in a ]
		samples = [ samples[i] for i in a ]
#  	filter  maf
	   	maf = get_MAF(gt)
		if(maf > maf_cut):
#	filter fisher test
			# ftest = fisherTest(gt)
			# if ftest > fisher_cut:
			a = [j for j in range(len(allsamples)) if allsamples[j] not in samples]
			for i in a:
				allgt[i] = -2 
			fout.write(key + "\t".join(map(str,allgt)) + "\n")
			# else :
				# cnt_f3 = cnt_f3 + 1
				# continue
   		else:
			cnt_f2 = cnt_f2 + 1
   			continue
	else:
		cnt_f1 = cnt_f1 + 1
		continue
fout2 = open(outpstat,'w')
fout2.write("Number_of_total_snp" + "\t" + "filtered_out1" + \
	"\t" + "filtered_out2" + \
	# "\t"+ "filtered_out3" + \
	"\t" + "Number_of_final")
fout2.write(nsnp,"\t" ,cnt_f1,"\t" ,cnt_f2,"\t" \
	# ,cnt_f3, "\t" \
	# ,nsnp - cnt_f1 - cnt_f2 - cnt_f3)
,nsnp - cnt_f1 - cnt_f2)
fout.close()
fout2.close()


