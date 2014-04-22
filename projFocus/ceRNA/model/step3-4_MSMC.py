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
# '''
argv     = sys.argv[1:]
debug    = False
# debug    = True
sltype   = "maxfreq" ## maxfreq, and other type underwork 
alpha    = 0.8
usage    = 'python ' + sys.argv[0] + '\n\
-a <alpha value: propotion of sample be covered> \n\
-m <mutataion matrix> \n\
-g <genesample list file>  \n\
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
        sltype    = arg
    elif opt in ("-o"):
        output      = arg
# '''
'''
mutfile     = "kegRegs_Apr-18-2014.summary.driverRegs.list.uniq.mut.matrix" 
gslistfile  = "gslist_Mar-24-2014_CnvMethSomFree.10smapMore.deg_20140325.txt.10more.hasReg.list" 
keyGeneFileDir = "/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/sigMut/step4_MSMC/test/test_sigKeygeneDir" 
output      = "sigGeneset_testSigKeyGeneDir.txt"
log         = output + ".log"
maxCombs    = 18
alpha       = 0.8
sltype      = 'minsize'
debug       = False
'''
##-----------------params
pvalCut     = 0.01 
effSSCutLow = 2
maxCombs    = 1000000
gssCut      = 1
output      = output + "_" + str(round(alpha,2)) + "_" + sltype 
log         = output + ".log"
logH        = open(log, 'wt')

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
outH.write("targetGene\tselectedSize\tselectedSet\n")
fileA  = [f for f in os.listdir(keyGeneFileDir) if f.endswith(".txt") ]
for f in fileA:
    start_time  = time.time()
    keygenefile = keyGeneFileDir+"/" + f
    mutDict, mutDictInfo = prepareDataMSMC(gslistfile, mutfile, keygenefile,\
                                           pvalCut = 0.01)
    if not mutDict or not mutDictInfo:
        logH.write("######" + f + "\n")
        logH.write("##r2.pval fail to  %f or no mutated regulator\n" %pvalCut )
        logH.write("##time for rurn(sec)\t" + str(time.time() - start_time) + "\n")
        continue
        sys.exit()
               
    allSolD, templog   = findAllSol( mutDictInfo['mutGintSmpNum'], mutDict, alpha = alpha,\
                ssCutLow = effSSCutLow, csCutUp = maxCombs,  gssCut = gssCut, log = True)
    resGSet   = selectSets(allSolD, type = sltype) 

    if debug:
        print mutDictInfo['tgene'] + "\t" + str(len(mutDictInfo['mutGintSmpNum']))+ \
                "\t" + str(len(mutDictInfo['mutRegs'])) 
        print mutDictInfo['tgene'],resGSet

    outH.write(mutDictInfo['tgene'] + "\t" + str(len(resGSet)) +\
               "\t" + ";".join(resGSet)+"\n")
    outH.flush()
    ##---writing log file
    logH.write("###----" + f + "\n")
    logH.write("##target_gene\t" + mutDictInfo['tgene']+"\n"  + \
               "##total_mutated_gint_sample\t"+ str(len(mutDictInfo['mutGintSmpNum']))+"\n"+ \
               "##total_sample\t"+ str(len(mutDictInfo['mutGintSmpNum']))+"\n"+ \
               "##total_mutated_driver\t" + str(len(mutDictInfo['mutRegs'])) + "\n")
    logH.write("#selection_type\t" + sltype + "\n" +\
                "#totoal_combination\t" + str(templog['totalCombCnt']) + "\n" +\
                "#effect_sample_size\t" + str(templog['effss']) + "\n" +\
                # "#sample size region\t" + templog['']+ "\n" +\
                "#minnimium_cost_selected\t" + str(templog['minCost']) + "\n"+\
                "#time_for_run(sec)\t" + str(time.time() - start_time) +"\n")
    logH.flush()
outH.close()
logH.close()
