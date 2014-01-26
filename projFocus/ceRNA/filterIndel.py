#!/usr/bin/env python
#input: 1.<file:vcf after gatk,annovar,and mutable filtering process >
#	2.<file:potential false positive gene names, one by one line>
#	3.<string: output file name> 
#output: 1.<file: VCF after filtering> 
#desp: filtering each final vcfs 

#loading package
import os
import sys,getopt
import re
import collections

#get input
# #-------test_begin-----
# inp = 'input_test.vcf'
# out = 'out_test.tsv'
# outlog = out + ".log"
# #------test_end-------

usage = "Usage: filterIndel.py -i <input.vcf> -o <out.vcf> \nExample: filterIndel.py -i /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/indelCall/test/TCGA-A1-A0SD-01A-11D-A10Y-09.bam_dindel_ouput.variantCalls.VCF -o /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/indelCall/test/TCGA-A1-A0SD-01A-11D-A10Y-09.bam_dindel_ouput.variantCalls.VCF.filtered.VCF"
argv = sys.argv[1:]
try:
      opts,args = getopt.getopt(argv,"hi:o:f:")
except getopt.GetoptError:
      print "Error Command line parameter \n" + usage 

for opt,arg in opts:
      if opt == "-h":
	    print usage 
	    sys.exit()
      elif opt == "-i":
            inp = arg
      elif opt == "-o":
            out = arg
            outlog = out + ".log"

#setting parameter
#gtnorm_cut = '0/0'
#pltum_cut = 40
#plnorm_cut = 30
#clr_cut = 1
dp3_cut = 1
dp4_cut = 1
dp34_cut = 3
dp_cut = 10
#mq0_cut =  0.1
#dpnorm_cut = 1
#mafnorm_cut = 0.05
#maftum_cut = 0
filter_val ='PASS'

outf = open(out,'w')
flag_header = 1 
cnt_line = 0 
cnt_res = 0
with open(inp) as inpf: 
      for line in inpf.readlines():
            if re.match("^#",line):
                  outf.write(line)
            else:
		   ##1       1276973 .       G       GACAC   379     PASS    DP=32;NF=9;NR=0;NRS=11;NFS=0;HP=1       GT:GQ   1/1:30
		   ###CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO=Depth of coverage; number of reads covering variants in the forward strand;number of reads covering....in the reverse strand; number of reads ... variants site in forward; number of reads ...variants site in reverse;    FORMAT  SAMPLE
                  cnt_line = cnt_line + 1
                  [chrom, pos, identifier, ref, alt, qual, filter, info, format, sample1] = \
                             line.strip().split("\t")
                  if filter == filter_val :
                        # print str(cnt_line) + "\t line processing..." 
                        ### getting info
                        infoArray = info.split(";")
                        #print "\t".join(infoArray) + "\n"  
                        infoDict = {}
                        flag=1
                        for infosub in infoArray:
                             temp = infosub.split("=")
                        #[tag,val] = infosub.split("=")
                             if len(temp) == 2:
                                   infoDict[temp[0]] = temp[1]
                             elif len(temp) == 1:
                                   infoDict[temp[0]] = ''
                             
                        if ('DP' in infoDict.keys()) and int(infoDict['DP']) < dp_cut:
                             flag = 0
                             #print "fail DP"
                             
                        if ('NFS' in infoDict.keys()) and int(infoDict['NFS']) < dp3_cut:
                             flag = 0
                             #print "fail NFS"
			if ('NRS' in infoDict.keys()) and int(infoDict['NRS']) < dp4_cut:
                             flag = 0  
			     #print "fail NRS"     
			if ('NRS' in infoDict.keys()) and ('NFS' in infoDict.keys()) and int(infoDict['NRS']) + int(infoDict['NFS']) < dp34_cut:
			     flag = 0
			     #print "fail total NFS, NRS"
                        ###-------------preparing for output 
                        if flag == 1:
                              cnt_res = cnt_res + 1
			      outf.write(line)
                        else:
                              continue
                  else :
                       continue

outlogf = open(outlog,'w')
print "Total recoreds processed: " + str(cnt_line) 
print "Recoreds survived after filtering: " + str(cnt_res)
outlogf.write("Total recoreds processed: " + str(cnt_line) +"\n") 
outlogf.write("Recoreds survived after filtering: " + str(cnt_res) +"\n")

print "#-----END------"

outf.close()
outlogf.close()
