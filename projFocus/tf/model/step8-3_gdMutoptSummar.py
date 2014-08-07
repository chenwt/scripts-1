#!/usr/bin/env python
# J.HE
# summarized the output from greedy into
# 1. selected sample matrix, target-sample
# 2. all actual mutation running results
#     a. statistics
#     b. tf-mutsample file
# 3. all actual & permutation results

import os, sys, re, getopt
from collections import defaultdict
import pandas as pd
import numpy as np
from pandas import DataFrame


###test begin---
# dataDir = "/Users/jh3283/projFocus/07252014/resultTFGreedy/data/test"
# ptAct = "_act_"
# ptPerm = "_pmLable"
# ptMatrix = "mutSampleVector"
# output = "summary_" 
###test end---

argv = sys.argv[1:]
usage = 'python ' + sys.argv[0] + '\
        -d <dataDir>  \
        -a <ptAct>  \
        -p <ptPerm>  \
        -m <ptMatrix>  \
        -o <output file prefix in working dir, 4 outputs, >'

example = 'python ' + sys.argv[0] + '\
        -d /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/tfMut/test  \
        -a _act_  \
        -p _pmLable  \
        -m _mutSampleVector  \
        -o test_summary'

try:
    opts,args = getopt.getopt(argv,"hd:a:p:m:o:")
except getopt.GetoptError:
    print usage + "\n" + example 
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print usage + "\n" + example 
        sys.exit()
    elif opt in ("-d"):
        dataDir = arg
    elif opt in ("-a"):
        ptAct = arg
    elif opt in ("-p"):
        ptPerm = arg
    elif opt in ("-m"):
        ptMatrix = arg
    elif opt in ("-o"):
        output = arg
###global variable

newOrder = [1, 5, 0, 4, 2, 3, 6, 7, 9 ,10, 8]
CNT=0
fgHd = 0


def getFileNames(dataDir):  
    tgenes = []; fileAct = []; filePm = []; fileMat = []
    for f in os.listdir(dataDir):
        if re.findall(ptMatrix, f):
            fileMat.append(f)
            tgenes.append(f.split("_")[2])
            continue 
        if re.findall(ptAct, f):
            fileAct.append(f)
        elif re.findall(ptPerm, f):
            filePm.append(f)
    return tgenes, fileAct, filePm, fileMat

def summaryMutSample(fileMat):
    '''
    merge all selected sample files
    further modify to make it do a summary
    '''
    outfh1 = open(output + ptMatrix, 'w'); fgHd = 0
    for f in fileMat:
        with(open(dataDir + "/" + f)) as fh:
            header = fh.readline()
            if fgHd  == 0:
                outfh1.write(header)
                fgHd = 1
            outfh1.write(fh.readline())

def recFormat(xLt):
    xLt = xLt.strip().split("\t")
    xSts =  map(int, [xLt[h] for h in newOrder[:4]]) + \
            map(float, [xLt[h] for h in newOrder[4:]] )
    xTms = xLt[11]
    return xSts, xTms

def summaryOpt(fileAct, t = 'act'):
    '''
    merge all actual mutation results
    further develop to run descriptive statistics
    '''
    actSumDt = defaultdict(list); actOptSmpDt = defaultdict() 
    for f in fileAct:
        with(open(dataDir + "/" + f)) as fh:             
            header = fh.readline().strip().split("\t")
            global fgHd 
            if fgHd == 0:
                header = header[1:]
                global headerSts 
                headerSts= [header[h] for h in newOrder ]
                global headerTms 
                headerTms = header[11]
                fgHd = 1
            line = fh.readline()
            if t == 'act':
                tgene, info = line.strip().split("\t",1)
                tempSum, tempTM = recFormat(info)
                actSumDt[tgene] = tempSum
                actOptSmpDt[tgene] = tempTM
                continue 
                
            cnt = ( int(f.split("_")[1][-1]) - 1) * int(f.split("_")[-1]) + 1
            while line:
                tgene, info = line.strip().split("\t",1)
                tempSum, tempTM = recFormat(info)
                actSumDt[tgene + "_"+ str(cnt) ] = tempSum
                actOptSmpDt[tgene + "_" + str(cnt) ] = tempTM
                cnt = cnt + 1
                line = fh.readline()
     
    actSumDF = DataFrame(actSumDt).T; actSumDF.columns = headerSts
    # print "all genes", len(fileAct)
    return actSumDF, actOptSmpDt


def summaryAllPerm(filePm):
    '''
    summarized permutation by compare to actual 
    1. before greedy p value
    2. after greedy p value
    3. greedy compareing p value
    '''
    allSumDF = DataFrame(np.zeros([1,len(actSumDF.columns)]), columns=actSumDF.columns)      
    
    for tg_crt in tgs:
        ipmf  = [f for f in filePm if re.findall("_"+tg_crt+"_", f)  and f]
        if len(ipmf) == 0:
            continue
        global CNT
        CNT = CNT + 1
        tempSum, _ = summaryOpt(ipmf, t = 'perm')
        allSumDF = pd.concat([allSumDF,tempSum])
    allSumDF = allSumDF.iloc[1:,:]
    print "permutated genes", CNT
    return allSumDF

def writeOptTgSmp(output, actOptSmpDt):
    '''
    given the dictionary of tgene: tf-smp; tfsmp, output it 
    '''
    outh = open(output + "_actOptTgSmp", 'w')
    outh.write("tgene\ttf_smp\n")
    for k,v in actOptSmpDt.items():
        outh.write(k + "\t" + v + "\n")
    outh.close()

### main -----
 
tgs, fileAct, filePm, fileSmp = getFileNames(dataDir)
actSumDF, actOptSmpDt = summaryOpt(fileAct)
permSumDF = summaryAllPerm(filePm)


### output data----
summaryMutSample(fileSmp)

actSumDF.to_csv(output + "_actSummary", sep= "\t")

writeOptTgSmp(output, actOptSmpDt)

permSumDF.to_csv(output + "_permSummary", sep= "\t")


### visulization for further develop 
def pyViz():
  import matplotlib.pyplot as plt
  
  actdata = actSumDF.loc[[tg_crt], ['pval_act','pval_opt', 'pval_actVopt']]
  ## actual compared to random
  
  permdata = permSumDF.loc[[g for g in permSumDF.index if re.findall(tg_crt, g)], \
                            ['pval_act','pval_opt','pval_actVopt']]
  permdata = -1 * log10(permdata)
  plt.figure(); 
  
  plt.scatter( permdata.iloc[:,1],  permdata.iloc[:,0])
  
  plt.figure(); 
  
  plt.scatter( permdata.iloc[:,1],  permdata.iloc[:,2])
  
print '[END]'
