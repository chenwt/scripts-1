#!/usr/bin/env python
#J.HE
#input: 1. <str: dirname including all the files>
#	2. <str: maf/gt>
#output:2. <file: chr_pos by patient-sample table with required information

import os
import glob
import sys
import re
import getopt
import string

argv = sys.argv[1:]
try:
  opts,args = getopt.getopt(argv,"hd:t:o:")
except getopt.GetoptError:
  print 'python genSomaticMutationForAllPatients.py -d <directory> -t <type maf/gt/raw/rawcount> -o <output file name> '
  sys.exit(2)

for opt, arg in opts:
  if opt == '-h':
    print 'python scripts_name.py -d <directory> -t <type maf/gt> -o <output file name> '
    sys.exit()
  elif opt in ("-d"):
    inpd      = arg
  elif opt in ("-t"):
    inptype   = arg
  elif opt in ("-o"):
    outp      = arg								

if inpd == '' or inptype == '' or outp == '' or not inptype in ['maf','gt','raw','rawcount']:
  print "input parameter ERROR!"
  print 'python filterSNP_utest_KWtest_v2.py -s <snp annotation file> -e <exp annotation file> -o <output file name> -j <number:cutoff for pval adjust 1e-6>'
  sys.exit(2)

print "working dir:\t" + inpd
print "information type:\t" + inptype
print "output file name:\t" + outp

#get all file name
fnames = ([file for root, dirs, files in os.walk(inpd)
      for file in files
          if file.endswith('.vcf') 
	  ])
outp  = os.getcwd() + "/" + outp
os.chdir(inpd)

#outpf	    = open(outp,'w') 
#
#if inptype == 'raw' or inptype == 'rawcount':
#    flag_header = 1 
#cnt_file    = 0
#
#for fname in fnames:
#  #pid.append(fname.split("."))
#  pid	    = fname.split(".")[0][:9]
#  cnt_file  = cnt_file + 1
#  chr_pos     = []
#  info_pid  = []
#  with open(fname) as f:
#    for line in f.readlines():
#      if not re.match(r"^#",line):
#	line_temp = line.strip().split("\t") 
#	#print "length of line:\t" + str(len(line_temp))
#	if inptype == "maf":
#	    [chrom, pos, info]  = [line_temp[0], line_temp[1], line_temp[7]]
#	    info	      = info.split(";")
#	    for a in info:  
#	      if re.match(r"^DP4",a):
#	          dp4	      = a
#		  dp4	      = map(int,string.replace(dp4,"DP4=","").split(","))
#		  info          = str(dp4[2] + dp4[3]) + "/" + str(dp4[0] + dp4[1])
#		#print "info:\t" + info
#		  info_pid.append(info)
#		  break	
#	elif inptype == "gt":
#	    [chrom, pos, info]  = [line_temp[0], line_temp[1], line_temp[9]]
#	    info	        = info.split(":")[0]
#	    info_pid.append(info)
#	elif inptype == 'raw':
#	    [chrom, pos, info]  = [line_temp[0], line_temp[1], line_temp[7]]
#	    for element in info.split(";"):
#	        if re.match(r"^I16",element):
#	            dp4  = map(int,string.replace(element,"I16=","").split(",")[:4])
#	            info = str(dp4[2] + dp4[3]) + "/" + str(dp4[0] + dp4[1])
#	            info_pid.append(info)
#	            break
#	#print pid + "\t" + info 
#	elif inptype == 'rawcount':
#	    [chrom, pos, info1, info2, info3] = [line_temp[0], line_temp[1],\
#                                                 line_temp[2],line_temp[3],\
#	                                         max(map(int,line_temp[4:]))]
#	    info = str(info3) + "/"+ str(int(info1) - int(info2))
#	    info_pid.append(info)
#	    
#	chr_pos.append(chrom + "_" + pos)
#  print "processsing file \t" + str(cnt_file) 
# 
#  if inptype == 'raw' or inptype == 'rawcount': 
#     if flag_header == 1:
#         outpf.write("patient_ID" + "\t" + "\t".join(chr_pos) + "\n")
#         flag_header = 0
#  else:
#      outpf.write("patient_ID" + "\t" + "\t".join(chr_pos) + "\n")    
#  outpf.write(pid + "\t" + "\t".join(info_pid) + "\n") 
#outpf.close()

#TODO
#post processing
if inptype == 'maf' or inptype == 'gt':
    chrPosAll = []
    #get all chr_pos number
    with open(outp) as f:
        for line in f.readlines():
            if re.match(r"^patient_ID",line):
                chr_pos = line.strip().split("\t")
                chrPosAll = list(set(chrPosAll + chr_pos))
                #numTemp = len(chr_pos)
                #if numTemp > numField:
                #    chrPosAll = chr_pos
                #    numField   = numTemp
    numField   = len(chrPosAll)
    #print numField
    #print chrPosAll
    #output data
    chrPosAll.sort()
    outputSort = open(outp + '.sort','w')
    outputSort.write("\t".join(chrPosAll) + "\n")
    
    with open(outp) as f:
        for i, line in enumerate(f):
            if re.match(r"^patient_ID",line):
                lineTemp = line.strip().split("\t")
                #print lineTemp
                nonZeroChrPosTemp = []
                for a in lineTemp :
                    if a in chrPosAll:
                        nonZeroChrPosTemp.append(chrPosAll.index(a))
                                
                #print nonZeroChrPosTemp
                #print len(nonZeroChrPosTemp)
                #if i == 9:
                #    sys.exit()
                next
            elif len(nonZeroChrPosTemp) > 0:
                info = ['NA'] * numField
                lineTemp = line.strip().split("\t")
                #print len(lineTemp)
                for (i,item) in zip(nonZeroChrPosTemp,lineTemp):
                    info[i] = item
                    #print i, item
                #print "\t".join(info)
                outputSort.write("\t".join(info) + "\n")
#                
    outputSort.close()
print "#-----END-------"
