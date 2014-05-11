#!/usr/bin/python
'''
Desp.: get summary for the result file of step3-4, 
'''

import os,sys
from collections import Counter, defaultdict
arg     = sys.argv[1:]
input   = arg[0]
outSum     = input + ".summary" 
outFreq     = input + ".summary.tarCnt" 

with(open(input)) as f:
    line = f.readline()
    line = f.readline()
    regTarD = defaultdict(list) 
    uniqReg = []
    while line:
        t, num, rs = line.strip().split("\t")
        regs    = rs.split(";")
        for r in regs:
            regTarD[r].append(t)
        uniqReg.extend(regs)
        line = f.readline()

regCount = Counter(uniqReg)
# topcut   = 10
# cutCount = 5
i        = 0
outFreqH = open(outFreq, 'w')
outFreqH.write("driverMutation\tnumTarget\ttargts\n")
for v, g in sorted(zip(regCount.values(),regCount.keys()), reverse = True):
    # if v > cutCount and i <= topcut :
    outFreqH.write(g + "\t" + str(len(regTarD[g])) + "\t"\
                       +  ":".join(regTarD[g]) + "\n" )
    
