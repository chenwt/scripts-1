#!/usr/bin/python
#J.HE

from    parseKeyRegFile import *
import  os, os.path 

# grpLasoResDir="/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/candiReg/run-Apr-1-2014/data"
grpLasoResDir="/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/candiReg/runApr30/data"
pvalCut = 0.01
output = "/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/candiReg/runApr30/summary/tgeneKeyreg" 
output  = output + "_" + str(pvalCut) 
outputH = open(output, 'w')
outputH.write("targetGene\tselectedCeRNADriver\n")
fileA = [f for f in os.listdir(grpLasoResDir) if f.endswith(".txt")] 
for f in fileA:
    file = grpLasoResDir +"/" + f
    print f
    t, crtRegs = parseKeyRegFile(file, pvalCut)
    if t and crtRegs:
        outputH.write(t[0] + "\t" + ";".join(crtRegs) + "\n")

outputH.close()

