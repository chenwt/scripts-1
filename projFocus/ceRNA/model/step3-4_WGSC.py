#!/usr/bin/python
#J.HE

'''
Desp.: weighted version of set covering.
output: selected set as well as covered samples
'''
import sys,getopt
import os, os.path
from parseKeyRegFile import parseKeyRegFile   
from parseGslistFile import parseGslistFile  
from prepareData_WGSC import *  
from wgsc import * 
import numpy as n
import time
from flattenSeq import flatten

# '''
argv     = sys.argv[1:]
debug    = False
alpha    = 0.85
wtype    = 'mean'
usage    = 'python ' + sys.argv[0] + '\n\
-a <alpha value: propotion of sample be covered> \n\
-m <mutataion matrix> \n\
-g <genesample list file>  \n\
-t <type for using zscore as weight>  \n\
-d <directory for keyRegulator result file> \n\
-o <output>'
try:
    opts,args = getopt.getopt(argv,"ha:m:g:d:t:o:")
except getopt.GetoptError:
    print usage 
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print usage  
        sys.exit()
    elif opt in ("-a"):
        alpha       = float(arg)
    elif opt in ("-m"):
        mutfile     = arg
    elif opt in ("-g"):
        gslistfile  = arg
    elif opt in ("-d"):
        keyGeneFileDir = arg
    elif opt in ("-t"):
        wtype    = arg
    elif opt in ("-o"):
        output      = arg
# '''

'''
mutfile     = "kegRegs_Apr-18-2014.summary.driverRegs.list.uniq.mut.matrix" 
gslistfile  = "gslist_Mar-24-2014_CnvMethSomFree.10smapMore.deg_20140325.txt.10more.hasReg.list" 
keyGeneFileDir = "/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/sigMut/step4_MSMC/test/test_sigKeygeneDir" 
output      = "sigGeneset_testSigKeyGeneDir.txt"
log         = output + ".log"
alpha       = 0.85
wtype      = 'mean'
debug       = True
'''
##-----------------params
pvalCut     = 0.01 
output      = output + "_" + str(round(alpha,2)) + "_" + wtype 
log         = output + ".log"
logH        = open(log, 'wt')
expTumorM   ="/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_tumor_Mar-21-2014.matrix_Mar-26-2014.voomNormed.matrix"
expNormalM  ="/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_normal_Mar-21-2014.matrix_Mar-26-2014.voomNormed.matrix"

def format2Dlist(l2d):
    out = [] 
    for i in l2d:
        out.append("["+",".join(map(str,i)) + "]")
    return ";".join(out)


logH.write('##Script path:\t'+ sys.argv[0]+"\n"+\
    "##mutation file\t" + mutfile +"\n"+\
    "##gene sample list file\t" + gslistfile + "\n"+\
    "##file directory " + keyGeneFileDir +"\n"+\
    "##output file" +  output +"\n"+\
    "##log_file " +  log + "\n"+\
    "##r2.pval_Cutoff\t" + str(pvalCut) + "\n"\
     "propotion Cutoff(alpha):\t" + str(alpha) + "\n")
     
##-----------------run set covering
outH   = open(output, 'w') 
logH.write("grpResFile" + "\t"+ "target_gene\t" + \
               "total_mutated_gint_sample\t" + \
               "total_sample\t" + \
               "total_mutated_driver\t" +\
               "selection_type\t"  +\
               "totoal_combination\t" +\
               "effect_sample_size\t" + \
               "minnimium_cost_selected\t" + \
               "time_for_run.sec\t" + "\n")

outH.write("tarGene\tnumIntSmp\tnumMutIntSmp\tnumCoverSmp\tnumReg\t\
           numDrReg\tnumMutDrReg\tnumIntMutDrReg\tnumSelectedReg\t\
           selectedGene\tselectedMutSmp\tselectedCost\n")
expnD  = loadNormExpfile(expNormalM)
zscoreD  = loadExpfile(expTumorM, expnD)

fileA  = [f for f in os.listdir(keyGeneFileDir) if f.endswith(".txt") ]

for f in fileA:
    start_time  = time.time()
    keygenefile = keyGeneFileDir+"/" + f
    crtTarRegMutObjD, mutDictInfo = prepareDataWGSC(mutfile, gslistfile, \
                                        keygenefile, zscoreD, pvalCut = 0.01)
    if not mutDictInfo:
        logH.write("##" + f + "\n")
        logH.write("##r2.pval fail to  %f or no mutated regulator\n" %pvalCut )
        logH.write("##time for rurn(sec)\t" + str(time.time() - start_time) + "\n")
        continue
    targ = crtTarRegMutObjD.keys()[0]
    if not crtTarRegMutObjD[targ][0]:
        continue
    geneL, smpL, costL   = wgsc( crtTarRegMutObjD[targ][0], \
                                crtTarRegMutObjD[targ][1], wtype = wtype, alpha = alpha)
    outNumberInfo = str(len(mutDictInfo['gintSmp'])) + "\t" +\
            str(len(mutDictInfo['mutGintSmp'])) + "\t" +\
            str(len(set(flatten(smpL)))) + "\t" + \
            str(len(mutDictInfo['allRegs'])) + "\t" + \
            str(len(mutDictInfo['mutRegs'])) + "\t" + \
            str(len(mutDictInfo['intMutRegs']))  + "\t" + \
            str(len(geneL)) 
    outH.write(targ + "\t" + outNumberInfo + "\t" + ";".join(geneL) + "\t" +\
               format2Dlist(smpL) + "\t" + ";".join(map(str,costL))+ "\n")
    outH.flush()

    logH.write(mutDictInfo['tgene'] + \
               str(time.time() - start_time) +"\n")
    logH.flush()
outH.close()
logH.close()
