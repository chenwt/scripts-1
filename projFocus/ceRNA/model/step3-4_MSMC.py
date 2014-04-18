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
##-----------------params
numIter     = 1000
result      = [ [] for _ in range(numIter)]
pvalCut     = 0.01 
alpha       = 0.8
outResFile  = output + '_' + str(alpha) + "_" + str(numIter)
outReslog   = outResFile + ".log"
logH        = open(outReslog, 'wt')

# print('Script path:\t'+ sys.argv[0])
# print('#Input files\ngene sample list file\t' + gslistfile +\
#       "\nkeyRegulatorfile\t" + keygenefile +\
#       "\nmutation matrix file\t" + mutfile)
# print('Output file:\t'+ output)
# print("#parameters:\npvalue for group lasso model\t %f\
#       \nportion of mutation gint sample\t %f\n " % (pvalCut, alpha))


##-----------------run MSMC
mutDict, mutDictInfo = prepareDataMSMC(gslistfile, mutfile, keygenefile,\
                                       pvalCut = 0.01)

# logH.write(outInfo['gene'] + "\t" + str(len(mutDictInfo['mutGintSmpNum']))+ \
        # "\t" + str(len(mutDictInfo['mutRegs'])) + "\n")
if not mutDict or not mutDictInfo:
    print "#not significant"
    # outputH = open(outResFile, 'w')
    # outputH.write("group lasso model r2.pval %s fail to pvalCutoff %f\n" % (mutDict, pvalCut) )
    # outputH.close()
    # logH.write("group lasso model r2.pval %s fail to pvalCutoff %f\n" % (mutDict, pvalCut) )
    sys.exit()

print mutDictInfo['tgene'] + "\t" + str(len(mutDictInfo['mutGintSmpNum']))+ \
        "\t" + str(len(mutDictInfo['mutRegs'])) 

sys.exit()
k           = int(alpha * len(mutDictInfo['mutGintSmpNum'] ) ) 
print "number of mutation gint samples:\t %d" % k

##----subsampling 
for i in range(numIter):
    '''# tempsample    = weighted_sample(mutdictinfo['mutgintsmpnum'], \
                mutdictinfo['mutgintsmpweight'],k) 
        ## weighted sampling methods
    '''
    tempsample     = random.sample(mutDictInfo['mutGintSmpNum'],k) 
    tempmutdict    = updateSets(mutDict, tempsample)
    tempressets, _ = msmc(tempmutdict.values(), set(tempsample))
    # print tempmutdict.keys(), set(tempsample)
    # print tempressets
    for gene, v in tempmutdict.items() :
        if v in tempressets:
            result[i].append(gene)   
    print "iteration " + str(i) + " selected genes " +  ";".join(result[i])

#----optimizing solution

out = defaultdict(int)
while result:
    gset = result.pop()
    g = ";".join(gset)
    if out.get(g,''):
        out[g][1] = out[g][1] + 1
    else:
        out[g] = [len(gset),1]

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

outputH = open(output + '_' + str(alpha) + "_" + str(numIter), 'w')
outputH.write(keyMax.replace(";","\n") + "\n")

outputH.close()
logH.close()


