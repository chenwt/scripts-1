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
from msmc import msmc, findMin, updateSets
from flattenSeq import flatten

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
print('Script path:\t'+ sys.argv[0])
print('Input file:\t' + input)
print('Output file:\t'+ output)

##-----------------run MSMC
mutDict, mutDictInfo = prepareDataMSMC(gslistfile, mutfile, keygenefile )
# print mutDict.items()
numIter = 1000
result  = [ [] for _ in range(numIter)]
alpha = 0.8
k     = int(alpha * len(mutDictInfo['mutGintSmpNum'] ) ) 
# print k
for i in range(numIter):
    # tempsample    = weighted_sample(mutdictinfo['mutgintsmpnum'], \
                    # mutdictinfo['mutgintsmpweight'],k) 
    tempsample     = random.sample(mutDictInfo['mutGintSmpNum'],k) 
    tempmutdict    = updateSets(mutDict, tempsample)
    # print tempmutdict.keys(), set(tempsample)
    tempressets, _ = msmc(tempmutdict.values(), set(tempsample))
    # print tempressets
    for gene, v in tempmutdict.items() :
        if v in tempressets:
            result[i].append(gene)   
    print "iteration " + str(i) + " selected genes " +  ";".join(result[i])

# out = Counter(result)
out = defaultdict(int)
while result:
    # for g in result.pop():
    g = ";".join(result.pop())
    if out.get(g,''):
        out[g] = out[g] + 1
    else:
        out[g] = 1

# if debug:
    # print sorted(out,key= out.get,reverse=True) 
print sorted(out.items(),key=lambda x:x[1],reverse=True) 

valMax = 0 
keyMax = []
for k, v in out.items():
    if v > valMax :
        valMax = v
        keyMax  = k  
    elif v == valMax and len(keyMax) > len(k):
        valMax = v
        keyMax  = k  
    else:
        continue

print keyMax, valMax

##-----------------prepare output 

# '''gene name, sample name, '''
outputH = open(output + '_' + str(alpha) + "_" + str(numIter), 'w')
# outputH.write("gene\tsamples\n")
outputH.write(keyMax.replace(";","\n") + "\n")

# # resGene = [ k for k, v in mutDict.items() if v in resSets]
# for k, v in mutDict.items() :
#     if v in resSets:
#         # print k, map(mutDictInfo['gintSamples'].__getitem__, map(lambda x:x-1, v))
#         outputH.write( k + "\t" + \
#                       ";".join(map(mutDictInfo['gintSamples'].__getitem__, \
#                                     map(lambda x:x-1, v))) + "\n" )

# outputH.close()

