#!/usr/bin/python
#J.HE
#input: <file: mutation file > <file: expression file>
#output: <file: mutation file for differential expression mutated genes>

import os
import sys, getopt
import math
import re
from scipy import stats
import numpy as np


usage = 'python step3-3_getFuncMutDriver.py -m <mutated gene file> -e <exp file\
> -o <output file name> -c <number:cutoff for pval adjust 1e-6>' 
error = "ERROR:"
argv = sys.argv[1:]
inps      = ''
inpe      = ''
pvalCut	  = 0.1  
outp      = ''
try:
    opts,args = getopt.getopt(argv,"hm:e:o:c:z:",\
                    ["mutfile=","expfile","ofile=","pval=", "debug="])
except getopt.GetoptError:
    print usage
    sys.exit(2)

debug = '' 
for opt, arg in opts:
    if opt == '-h':
        print usage 
        sys.exit()
    elif opt in ("-m","--mutfile"):
        mutfile = arg
    elif opt in ("-e","--expfile"):
        expfile = arg
    elif opt in ("-c","--pval"):
        pvalCut = float(arg)
    elif opt in ("-o","--ofile"):
        output = arg
    elif opt in ("-z","--debug"):
        debug = arg

outputMut = output + "_" + str(pvalCut) + ".matrix" 
outputStat = output + "_" + str(pvalCut) + ".stats" 
outlog  = output + "_" + str(pvalCut)  + ".log"
outlogf = open(outlog,'w')


if mutfile == '' or expfile == '' or output == '':
    print usage
    sys.exit(2)

print('Input file:\t' + mutfile)
print('Database file:\t' + expfile)
print('outputut file:\t'+ output + ".adjPass_" + str(pvalCut) + ".mat" + "\t" + output + ".adj.snp")
print('debug:\t'+ debug )

outlogf.write('Input file:\t'     + mutfile	  + "\n")
outlogf.write('Database file:\t'  + expfile	  + "\n")
outlogf.write('outputut file:\t'  + output	  + "\n")
# outlogf.write('Log file:\t'       + outlog    + "\n")

##----------------------------
def rUTest(x,y):
    if len(x) != len(y):
      print "ERROR: different lenght for Utest"
      sys.exit(2)
    else:
      z , p = stats.wilcoxon(x, y)
      # p = r.wilcox_test(x,y)['p.value']
      # z = r.wilcox_test(x,y)['statistic']['W']
    return([z,'%e' % float(p)])
def tTest(x,y):
    return stats.mstats.ttest_ind(x, y) 
def kwTest(x,y):
    [t,p] = stats.mstats.kruskalwallis(x,y)
    return([t,'%e' % float(p)])
def myZscore(x,y):
    zs = map(lambda xx:( xx - np.mean(y)) / np.std(y), x)
    z  = sum(zs) / np.sqrt(len(zs))
    p  = stats.norm.sf(z) 
    return  z, p 

def pAdj(plist,type):
     plist = map(float,plist)
     plist_sorted = sorted(plist)
     m = len(plist)
     p_new = []
     if   type == 'b':
       for i, pval in enumerate(plist):
           p_new.append(min(pval * m / (plist_sorted.index(pval) + 1),1))
     elif type == 'f':
       """ need to develop """
       pass
     return(p_new) 

def z2p(z):
    from math import erf,sqrt
    p=0.5*(1+erf(z/sqrt(2)))
    return (1-p)

def getStats(mutArray, expArray, smpmut, smpexp) :
    idExp = np.in1d(smpexp, smpmut[ np.nonzero(mutArray > 0)]) 
    if len(idExp) >0 : 
        idnotExp = np.in1d(smpexp, smpmut[ np.nonzero(mutArray == 0)])
        return myZscore(expArray[idExp], expArray[idnotExp])
    else:
        return 0.0, 1

mutDict = {}
#----load mutation data
with(open(mutfile)) as f:
    smpMutArray = np.array(map(lambda x:x[5:16].replace("-","."),\
                              f.readline().strip().split("\t")[1:] ))
    line = f.readline()
    while line:
        gene, val = line.strip().split("\t",1)
        mutDict[gene] = np.array(map(int,val.split("\t")))
        line = f.readline()
outlogf.write("number of input mutation:\t" + str(len(mutDict.keys())) + "\n")
outputMutH  = open(outputMut, 'wt') 
outputMutH.write("gene" + "\t" + "\t".join(smpMutArray.tolist()) + "\n")
outputStatH = open(outputStat, 'wt') 
outputStatH.write("gene\tzscore\tpval\n")

cnt = 0 
with open(expfile) as f:
    line = f.readline()
    smpExpArray = np.array(line.strip().split("\t"))
    line = f.readline()
    while line:
        gene, expval   = line.strip().split("\t",1)
        expval         = np.array(map(float, expval.split("\t")))
        if len(mutDict.get(gene,'')) > 0:
            if debug:
                mutVal = mutDict.get(gene,'') 
            zscore, pval = getStats(mutDict[gene], expval,\
                                    smpMutArray, smpExpArray ) 
            if pval <= pvalCut:
                cnt = cnt  + 1
                outputMutH.write(gene + "\t" + "\t".join(map(str, \
                    mutDict[gene].tolist())) + "\n")
                outputStatH.write("\t".join(map(str,[gene, zscore, pval])) + "\n")
                if debug:
                    print gene, zscore, pval 
        line = f.readline()

outlogf.write("number of functional mutations:\t" + str(cnt) + "\n")

outlogf.close()
outputMutH.close()
outputStatH.close()
