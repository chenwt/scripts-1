#!/usr/bin/python
#J.HE
'''
Desp.: given result from candiReg, gslist, mutations matrix, find the minium
size of driver regulator set which are mutated in most of Gint samples 
imput: matrix of mutation
output: selected set as well as covered samples
'''

import sys,getopt
import re
from collections import defaultdict
from parseKeyRegFile import parseKeyRegFile 
from parseGslistFile import parseGslistFile  
import numpy as np
from   collections import Counter, Sequence
from flattenSeq import flatten
# argv    = sys.argv[1:]
# input   = ''
# output  = ''
# usage = 'python ' + sys.argv[0] + ' -i <input>  -o <output>'
# example = 'python ' + sys.argv[0] + ' -i <input>  -o <output>'
# try:
#     opts,args = getopt.getopt(argv,"hm:g:c:o:")
# except getopt.GetoptError:
#     print usage + "\n" + example 
#     sys.exit(2)
# for opt, arg in opts:
#     if opt == '-h':
#         print usage + "\n" + example 
#         sys.exit()
#     elif opt in ("-m"):
#         mutfile = arg
#     elif opt in ("-g"):
#         gslistfile = arg
#     elif opt in ("-c"):
#         keygenefile = arg
#     elif opt in ("-o"):
#         output = arg
# print('Script path:\t'+ sys.argv[0])
# print('Input file:\t' + mutfile + "\t" + gslistfile + "\t" + keygenefile)
# print('Output file:\t'+ output)

def formatSampleName(code19):
    if len(code19) >11:
        return code19[5:16].replace("-",".")
    else :
        return code19.replace("-", ".")

def unionDictList(d):
    out = set() 
    for k, v in d.items():
       out.update(v) 
    return out

keyRegMutDict = defaultdict(set)
def prepareDataMSMC( gslistfile, mutfile, keygenefile, pvalCut = 0.01, debug = False ):
    tgeneSum, regsSum   = parseKeyRegFile(keygenefile, pvalCut)
    reglist             = regsSum
    tgene = tgeneSum[0]
    gintsmplist         = parseGslistFile(tgene, gslistfile)
    if debug:
        print "Debug-------"
        print tgene
        print reglist
        print len(gintsmplist), gintsmplist[:5]
        print len(reglist), reglist[:5]
        print "Debug-------"
    with(open(mutfile)) as f:
         gene, allsamples = f.readline().strip().split("\t",1) 
         allsamples       = allsamples.split("\t")  
         allsamples       = map(formatSampleName, allsamples)
         gIntIndex        = [id for id, a in enumerate(allsamples) if a in gintsmplist]
         gIntSmps         = map(allsamples.__getitem__, gIntIndex)
         if debug:
             print len(allsamples), allsamples[:5]
             print len(allsamples), allsamples[:5]
             print len(gIntIndex),  gIntIndex[:5]
         line = f.readline()
         while line:
             gene, vals      = line.strip().split("\t", 1)
             vals = map(int, vals.split("\t"))
             if gene in reglist and reduce(lambda x,y: x + y,\
                                           map(vals.__getitem__, gIntIndex )) > 0 :
                    keyRegMutDict[gene] = map(vals.__getitem__, gIntIndex )      
                    keyRegMutDict[gene] = set([idx + 1 for (idx, v) in \
                                enumerate(keyRegMutDict[gene]) if v != 0 ])
                    # print gene,keyRegMutDict[gene], map(gIntSmps.__getitem__, map(lambda x:x-1, keyRegMutDict[gene]))  
             line        = f.readline()
    outInfo = {}
    tempCounter = Counter(list(flatten(keyRegMutDict.values())))
    outInfo['numGintSmp']       = len(gIntIndex)
    outInfo['gintSamples']      = gIntSmps 
    outInfo['mutGintSmpNum']    = tempCounter.keys()
    outInfo['mutGintSmpWeight'] = tempCounter.values()
    outInfo['tgene']            = tgene
    outInfo['mutGintSmp']       = map(gIntSmps.__getitem__,\
                                   map(lambda x:x-1, outInfo['mutGintSmpNum']))  
    return  keyRegMutDict, outInfo  
def __test__():
    prepareDataMSMC(gslistfile, mutfile, keygenefile, debug = True)

def __main__():
    print "Run test"
    __test__()

