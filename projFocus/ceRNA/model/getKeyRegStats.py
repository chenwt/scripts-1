#!/usr/bin/python
#J.HE
#Desp.: Given the result for getKeyReg.r, gentrate the data for presentation
#input: file direct or/ dirs
#output:

import os, re
import sys,getopt
import numpy as np
import textwrap
argv = sys.argv[1:]
input = ''
output = ''
usage = 'python getKeyRegStats.py -d <input file directory> \
        -o <output file name>'
example = 'python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/getKeyRegStats.py\
-d /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/candiReg/test\
-o test.out'
try:
    opts,args = getopt.getopt(argv,"hd:o:")
except getopt.GetoptError:
    print(textwrap.fill(usage + "\n" + example, 60))
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print usage + "\n" + example 
        sys.exit()
    elif opt in ("-d"):
        fileDir = arg
    elif opt in ("-o"):
        output = arg
print('Script path:\t'+ sys.argv[0])
print('Input file:\t' + input)
print('Output file:\t'+ output)

outputSig = output + ".sigcnt"
pval_cut = 0.01 
def parseFile(file, pval_cut):
    print file
    with(open(file)) as f:
        line = f.readline()
        tgene = file.strip().split("_",1)[0]
        tsum = [tgene, 0.0, 0.0, 0, 0]
        for i, line in enumerate(f):
            if i >= 5:
                break
            else:
                if re.match(r"^#r2\t", line):
                   tsum[1] = line.strip().split("\t")[1] 
                elif re.findall(r"^#r2.pval", line):
                    _, pval = line.strip().split("\t")
                    if float(pval) <= pval_cut:
                        tsum[2] = str(pval)
                    else:
                        tsum = '' 
                        break
                elif re.match(r"^#totalReg", line):
                    tsum[3] = line.strip().split("\t")[1]
                elif re.match(r"^#sigReg",line):
                    tsum[4] = line.strip().split("\t")[1]
                elif re.match(r"^#target", line):
                    if not line.strip().split("\t")[1] == tgene:
                        print "error at gene name matching"
                else:
                    pass
    return tsum

fArray = os.listdir(fileDir)
outheader = "\t".join(['gene_tar','r2', 'r2pval', 'regTotal', 'regSig']) 
sufix = '.txt'
cntSig = 0 
outputH = open(output, 'w')
outputH.write(outheader + "\n")
for name in [name for name in fArray if name.endswith(sufix)]:
    out = parseFile(name, pval_cut)
    if out :
        cntSig = cntSig + 1
        outputH.write("\t".join(map(str,out)) + "\n")
open(outputSig,'w').write("totalTargets\t" + str(len(fArray)) + \
                          "\tsigTarget\t" + str(cntSig) + "\n")
outputH.close()

