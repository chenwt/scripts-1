#!/usr/bin/python
#J.HE
#Desp: replacement of shell script < grep -wf snpid.txt snp.mat.anno > since it is TOO SLOW!!
#	python only need several SECONDS!!!
#input: 1. <ID file: one each raw> 
#       2. < annotated ID mat file: first column is ID> 
#output: 2. < snp mat file of te given snp> 


import os,sys,getopt

argv = sys.argv[1:]
inputID  = ''
inputMat = ''
output   = ''

usage = '~/tools/python/Python_current/python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/extractSNPmat.py -i ~/SCRATCH/projFocus/ceRNA/knowledgeBase/GWAS_catalog_brca_SNPid.txt -m ../brca_snp_tumor_731.mat.anno -o brca_snp_tumor_731_GWASCatalogSNP.mat.anno'
example = 'python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/extractSNPmat.py -i <input snpid> -m <input snp annoataed matrix> -o <output>'
  
try:
  opts,args = getopt.getopt(argv,"hi:m:o:")
except getopt.GetoptError:
  print usage + "\n" + example
  sys.exit()

for opt, arg in opts:
    if opt == '-h':
          print usage +"\n" + example
	  sys.exit()
    elif opt in ("-i"):
	  inputID = arg
    elif opt in ("-m"):
          inputMat = arg
    elif opt in ("-o"):
          output = arg

# print inputID
# print inputMat
# print output
snpid = []

with open(inputID) as fhID:
  for line in fhID.readlines():
    snpid.append(line.strip())
# print "snpid number: \t" + str(len(snpid)) 
	
fhout  = open(output,'w') 
cntOut = 0 
with open(inputMat) as fhMat:
  for i, line in enumerate(fhMat):
    if i == 0:
        fhout.write(line)
    else :
	[key, val] = line.split("\t",1)
	for eleSnpid in snpid:
	  if eleSnpid == key:
	    # print eleSnpid
	    fhout.write(line)
	    cntOut = cntOut + 1
	  else :
	    continue
   
print "output snpid records: \t" + str(cntOut)
print "#-----END----"

