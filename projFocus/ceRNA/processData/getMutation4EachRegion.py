#!/usr/bin/env python
#input: 1. <file:file of matrix of mutations in difference samples> 
#       2. <file: gene, and location of interest:format: identifier: chr. pos. strand>
#output: 1. <file: .mat.anno file including all mutation info for the input files>

USAGE  ="Usage: python getMutation4EachRegion.py  \n\
        -c <file:full path of cnv level3 files,one each line> \n\
        -e <file: gene_exp.mat> \n\
        -r <overlap region size: 2000 bp>\n 
        -o <filename: output file name>"
EXAMPLE="Example: python getMutation4EachRegion.py  \n\
        -c  \n\
        -e  \n\
        -r 2000 \n 
        -o "

import os, re
import sys, getopt

argv = sys.argv[1:]
try:
    opts,args = getopt.getopt(argv, "h:c:e:o:")
except:
    print USAGE 
    sys.exit()
for opt,arg in opts:
    if opt == '-h':
      print USAGE 
      print EXAMPLE
      sys.exit()
    elif opt == '-m':
      inputmaf = arg
    elif opt == '-c':
      region_cut = int(arg)
    elif opt == '-r':
      inpRegion = arg
    elif opt == '-o':
      output = arg
      outputfreq = output + ".freq"
      outputlog = out + ".log"

outlogf = open(outlog,'w')

##-----setting parameters
##region cut for tss 2000, 3'UTR 0, 5'UTR 0 
# region_cut = 2000 

def chrom2Num(chrom):
    if chrom == 'X':
        chrom = 23
    if chrom == 'Y':
        chrom = 24
    chrom = int(chrom) - 1
    return(chrom)

def num2Chrom(num):
    num = int(num)
    if num == 22:
       num = 'X' 
    elif num == 23:
      num = 'Y' 
    else:
      num = num + 1
    return(str(num))

###----loading_Gene_Information----
allChroms = map(str,range(22)) + ['X','Y', 'x', 'y']
chromArray = [[] for _ in range(24)]
genePosDict  = {}
cnt = 0
outputfreqH = open(outputfreq, 'w+')

with open(inpRegion) as f:
    line = f.readline() 
    while line
        cnt = cnt + 1
        if re.match("gene",line):
            outputfreqH.write(line.strip().split("\t",3) \
                              + "mutationFreq" + "\n")
        else:
            [identifer, chrom, pS, pE, strand] = \
                    line.strip().split("\t")[:5]
            if chrom  in allChroms: 
	            genePosDict[(chrom,pS,pE)] = identifer
	            chrom = chrom2Num(chrom)
	            chromArray[chrom].append([int(pS), int(pE)])
        line = f.readline()

outlogf.write("number of input Region\t" + str(cnt) + "\n")
for chromPosArray in chromArray:
    chromPosArray.sort(key=lambda x: x[0])
# chromPosArray.sort(key=lambda tup:tup[0])
## print "chrom\t" + str(chromArray.index(chromPosArray)) +"\tgenes:\t" + str(len(chromPosArray)) 

#####--------------------------------loading maf. 
outputH = open(output, 'w+')

cnt = 0 
# chromDict = [{} for _ in range(24)]
with open(inputmaf) as f:
    line = f.readline()
    while line
        cnt = cnt + 1
        if re.match("^gene",line):
            outputH.write(line)
        else:
            [identifier,chrom,pStart,pEnd,val] = line.strip().split("\t",4)
            if chrom in allChroms:
                chrom = chrom2Num(chrom)
                tempStart = long(pStart) - region_cut
                tempEnd   = long(pEnd)   + region_cut
                for [pps,ppe] in [[pps, ppe] for [pps,ppe] in chromArray[chrom] \
                            if pps >= tempStart and ppe <= tempEnd]:
                    ### output
                    outRecord = genePosDict[(chrom,pps,ppe)] + "\t" + identifier\
                    + chrom + pStart + pEnd + val + "\n" 
                    print outRecord
                    print outfreqRecord
                    outputH.write(outRecord)
                    outputfreqH.write(outfreqRecord)
            else:
                continue
            outfreqRecord = genePosDict[(chrom, pps, ppe)] + \
                            "\t".join([chrom, pps, ppe)] + "\n"
            print outfreqRecord  
            outputfreqH.write(outfreqRecord)
        line = f.readline()
outlogf.write("total\t" + str(cnt) + " association founded\n")  

outputH.close()
outputfreqH.close()
print "#--------------DONE--------"
