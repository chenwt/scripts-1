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
  print 'python scripts_name.py -d <directory> -t <type maf/gt> -o <output file name> '
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

if inpd == '' or inptype == '' or outp == '':
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
os.chdir(inpd)

outpf	    = open(outp,'w') 
flag_header = 1 
cnt_file    = 0

for fname in fnames:
  #pid.append(fname.split("."))
  pid	    = fname.split(".")[0][:9]
  cnt_file  = cnt_file + 1
  chr_pos     = []
  info_pid  = []
  with open(fname) as f:
    for line in f.readlines():
      if not re.match(r"^#",line):
	line_temp = line.strip().split("\t") 
	#print "length of line:\t" + str(len(line_temp))
	if inptype == "maf":
	    [chr, pos, info]  = [line_temp[0], line_temp[1], line_temp[7]]
	    info	      = info.split(";")
	    for a in info:  
	      if re.match(r"^DP4",a):
		dp4	      = a
		dp4	      = map(int,string.replace(dp4,"DP4=","").split(","))
		info          = str(dp4[0] + dp4[1]) + "/" + str(dp4[2] + dp4[3])
		#print "info:\t" + info
		info_pid.append(info)
		break	
	elif inptype == "gt":
	    [chr, pos, info]  = [line_temp[0], line_temp[1], line_temp[9]]
	    info	      = info.split(":")[0]
	    info_pid.append(info)
	#print pid + "\t" + info 
	chr_pos.append(chr + "_" + pos)
  print "processsing file \t" + str(cnt_file) 
  outpf.write("patient_ID" + "\t" + "\t".join(chr_pos) + "\n")
  outpf.write(pid + "\t" + "\t".join(info_pid) + "\n") 
outpf.close()

#post processing

print "#-----END-------"
