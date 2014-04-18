#!/usr/bin/python
#J.HE
'''
Desp.: given result from candiReg, gslist, mutations matrix, find the minium
size of driver regulator set which are mutated in most of Gint samples 
imput: matrix of mutation
output: selected set as well as covered samples
'''

import sys,getopt
from collections import defaultdict, Counter
import re
import random
from parseKeyRegFile import parseKeyRegFile   
from parseGslistFile import parseGslistFile  
from prepareData_MSMC import prepareDataMSMC  
from msmc import * 
from flattenSeq import flatten

def getCmdArgs():
    argv     = sys.argv[1:]
    input    = ''
    output   = ''
    debug    = True
    usage    = 'python ' + sys.argv[0] + ' -m <mutataion matrix>  -g <genesample \
    list file -c <keyRegulator result file> -o <output>'
    example  = 'python ' + sys.argv[0] + ' -i <input>  -o <output>'
    try:
        opts,args = getopt.getopt(argv,"hm:g:c:o:")
    except getopt.GetoptError:
        print usage + "\n" + example 
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print usage + "\n" + example 
            sys.exit()
        elif opt in ("-m"):
            mutfile = arg
        elif opt in ("-g"):
            gslistfile = arg
        elif opt in ("-c"):
            keygenefile = arg
        elif opt in ("-o"):
            output = arg

mutfile     = "kegRegs_Apr-18-2014.summary.driverRegs.list.uniq.mut.matrix" 
gslistfile  = "gslist_Mar-24-2014_CnvMethSomFree.10smapMore.deg_20140325.txt.10more.hasReg.list" 
keygenefile = "RNF11_candidateRegs_Mar-31-2014.txt" 
output      = "RNF11.sigGeneset"
log         = output + ".log"

##-----------------params
numIter     = 1000
pvalCut     = 0.01 
alpha       = 0.8
logH        = open(log, 'wt')
debug       = True

if debug:
    print "mutation file", mutfile
    print "gene sample list file", gslistfile
    print "result from group lasso ", keygenefile
    print "output file", output
    print "log file ", log

if debug:
    print('Script path:\t'+ sys.argv[0])
    print('#Input files\ngene sample list file\t' + gslistfile +\
          "\nkeyRegulatorfile\t" + keygenefile +\
          "\nmutation matrix file\t" + mutfile)
    print('Output file:\t'+ output)
    print("#parameters:\npvalue for group lasso model\t %f\
          \nportion of mutation gint sample\t %f\n " % (pvalCut, alpha))
    
##-----------------run MSMC
mutDict, mutDictInfo = prepareDataMSMC(gslistfile, mutfile, keygenefile,\
                                       pvalCut = 0.01)

if not mutDict or not mutDictInfo:
    print "#not significant"
    # outputH = open(outResFile, 'w')
    # outputH.write("group lasso model r2.pval %s fail to pvalCutoff %f\n" % (mutDict, pvalCut) )
    # outputH.close()
    # logH.write("group lasso model r2.pval %s fail to pvalCutoff %f\n" % (mutDict, pvalCut) )
    sys.exit()

if debug:
    print mutDictInfo['tgene'] + "\t" + str(len(mutDictInfo['mutGintSmpNum']))+ \
            "\t" + str(len(mutDictInfo['mutRegs'])) 
    print mutDict.items()
    print mutDictInfo['mutGintSmpNum']
    # print "number of mutation gint samples:\t %d" % k
    # logH.write(mutDictInfo['tgene'] + "\t" + str(len(mutDictInfo['mutGintSmpNum']))+ \
            # "\t" + str(len(mutDictInfo['mutRegs'])) + "\n")

allSolD  = findAllSol( mutDictInfo['mutGintSmpNum'], mutDict, debug = True) 
resGSet  = selectSets(allSolD)

if debug:
    print "all solution\t", allSolD
    print "selected sets\t",resGSet
##-----------------prepare output 

logH.close()


