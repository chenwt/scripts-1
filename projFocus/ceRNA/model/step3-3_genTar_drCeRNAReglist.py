#!/usr/bin/python
#J.HE

from    parseKeyRegFile import *
import  os, os.path 

import sys,getopt
argv = sys.argv[1:]
input = ''
output = ''

try:
    opts,args = getopt.getopt(argv,"hi:o:")
except getopt.GetoptError:
    print usage + "\n" + example 
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print usage + "\n" + example 
        sys.exit()
    elif opt in ("-i"):
        grpLasoResDir = arg
    elif opt in ("-o"):
        output = arg
print('Script path:\t'+ sys.argv[0])
print('Input file:\t' + input)
print('Output file:\t'+ output)

# grpLasoResDir="/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/candiReg/run-Apr-1-2014/data"
# grpLasoResDir="/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/candiReg/runApr30/data"
pvalCut = 0.01
# output = "/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/candiReg/runApr30/summary/tgeneKeyreg" 
output  = output + "_" + str(pvalCut) 
outputH = open(output, 'w')
outputH.write("targetGene\tselectedCeRNADriver\n")
# fileA = [f for f in os.listdir(grpLasoResDir) if f.endswith(".txt")] 
fileA = [f for f in os.listdir(grpLasoResDir) if f.endswith(".txt.sigreg")] 
for f in fileA:
    file = grpLasoResDir +"/" + f
    print f
    t, crtRegs = parseKeyRegFile(file, pvalCut)
    if t and crtRegs:
        outputH.write(t[0] + "\t" + ";".join(crtRegs) + "\n")

outputH.close()

