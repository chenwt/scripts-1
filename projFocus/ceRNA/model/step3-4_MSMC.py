#!/usr/bin/python
#J.HE
'''
Desp.: given result from candiReg, gslist, mutations matrix, find the minium
size of driver regulator set which are mutated in most of Gint samples 
imput: matrix of mutation
output: selected set as well as covered samples
'''

import sys,getopt
import os, os.path
from parseKeyRegFile import parseKeyRegFile   
from parseGslistFile import parseGslistFile  
from prepareData_MSMC import prepareDataMSMC  
from msmc import * 
import time

argv     = sys.argv[1:]
debug    = False
usage    = 'python ' + sys.argv[0] + '\n\
-m <mutataion matrix> \n\
-g <genesample list file>  \n\
-d <directory for keyRegulator result file> \n\
-o <output>'
try:
    opts,args = getopt.getopt(argv,"hm:g:d:o:")
except getopt.GetoptError:
    print usage 
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print usage  
        sys.exit()
    elif opt in ("-m"):
        mutfile     = arg
    elif opt in ("-g"):
        gslistfile  = arg
    elif opt in ("-d"):
        keyGeneFileDir = arg
    elif opt in ("-o"):
        output      = arg
        log         = output + ".log"
'''
mutfile     = "kegRegs_Apr-18-2014.summary.driverRegs.list.uniq.mut.matrix" 
gslistfile  = "gslist_Mar-24-2014_CnvMethSomFree.10smapMore.deg_20140325.txt.10more.hasReg.list" 
keyGeneFileDir = "/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/sigMut/step4_MSMC/test/test_sigKeygeneDir" 
output      = "sigGeneset_testSigKeyGeneDir.txt"
log         = output + ".log"
debug       = True
'''
##-----------------params
pvalCut     = 0.01 
alpha       = 0.8
logH        = open(log, 'wt')
sltype      = "minsize" ## maxfreq, and other type underwork 
effSSCutLow = 5
maxCombs    = 1000000


if debug:
    print('Script path:\t'+ sys.argv[0])
    print "##mutation file", mutfile
    print "##gene sample list file", gslistfile
    print "##file directory ", keyGeneFileDir
    print "##output file", output
    print "##log file ", log

logH.write('##==Script path:\t'+ sys.argv[0]+"\n"+\
    "##==mutation file\t" + mutfile +"\n"+\
    "##==gene sample list file\t" + gslistfile + "\n"+\
    "##==file directory " + keyGeneFileDir +"\n"+\
    "##==output file" +  output +"\n"+\
    "##==log file " +  log + "\n"+\
    "##==parameters:\nr2.pval Cutoff\t %f\
          \npropotion Cutoff(alpha):\t %f\n\
        minEffectSampleSize\t%i\nmaxCombination:\t%i\n"\
            % (pvalCut, alpha, effSSCutLow, maxCombs))
    

if debug:
    # print('#Input files\ngene sample list file\t' + gslistfile +\
          # "\nkeyRegulatorfile\t" + keygenefile +\
          # "\nmutation matrix file\t" + mutfile)
    print('Output file:\t'+ output)
    print("#parameters:\npvalue for group lasso model\t %f\
          \nportion of mutation gint sample\t %f\n " % (pvalCut, alpha))
    
##-----------------run MSMC
outH   = open(output, 'w') 
outH.write("targetGene\tnumSelectedSet\tselectedSet\n")
fileA  = [f for f in os.listdir(keyGeneFileDir) if f.endswith(".txt") ]
for f in fileA:
    start_time  = time.time()
    keygenefile = keyGeneFileDir+"/" + f
    mutDict, mutDictInfo = prepareDataMSMC(gslistfile, mutfile, keygenefile,\
                                           pvalCut = 0.01)
    if not mutDict or not mutDictInfo:
        logH.write("######" + keygenefile + "\n")
        logH.write("##r2.pval fail to  %f or no mutated regulator\n" %pvalCut )
        logH.write("##time for rurn(sec)\t" + str(time.time() - start_time) + "\n")
        continue
        sys.exit()
    
               
    allSolD   = findAllSol( mutDictInfo['mutGintSmpNum'], mutDict, alpha = alpha,\
                ssCutLow = 5, csCutUp = 1000000, gssCut = 1)
    resGSet   = selectSets(allSolD, type = sltype) 

    if debug:
        print mutDictInfo['tgene'] + "\t" + str(len(mutDictInfo['mutGintSmpNum']))+ \
                "\t" + str(len(mutDictInfo['mutRegs'])) 
        print mutDictInfo['tgene'],resGSet

    outH.write(mutDictInfo['tgene'] + "\t" + str(len(resGSet)) +\
               "\t" + "||".join(resGSet)+"\n")
    outH.flush()
    ##---writing log file
    logH.write("######" + f + "\n")
    logH.write("##target gene\tmutDictInfo['tgene']\n"  + \
               "##total mutated gint sample\t"+ str(len(mutDictInfo['mutGintSmpNum']))+"\n"+ \
               "##total sample\t"+ str(len(mutDictInfo['mutGintSmpNum']))+"\n"+ \
               "##total mutated driver\t" + str(len(mutDictInfo['mutRegs'])) + "\n")
    
    logH.write( "#effect samplesize\t" +" "+ "\n" +\
                "#totoal combination\t" +" " + "\n" +\
                "#sample size region\t" + "" + "\n" +\
                "#selection type\t" + "" + "\n" +\
                "#minnimium cost selected\t" + "" + "\n"+\
                "#time for run(sec)\t" + str(time.time() - start_time) +"\n")
    logH.flush()
outH.close()
logH.close()
