import  random
import  bisect
from    flattenSeq  import flatten
from    collections import Counter, Sequence, defaultdict
from    itertools   import combinations
import  scipy.misc
import  math
import  itertools

def findMin(S, R):
    # S.sort(key=lambda x:len(x),reverse=True)
    minCost = 99999.0
    w = [1] * len(S)
    minElement = -1
    lenS = map(len, S)
    for i, s in enumerate(S):
        try:
            cost = w[i]/(len(s.intersection(R)))
            if cost < minCost:
                minCost = cost
                minElement = i
        except:
            pass
    return S[minElement], w[minElement]

def msmc(S, U, alpha = 0.8 ):
    R = U
    C = []
    costs = []
    while len(R) != 0:
       S_i, cost = findMin(S, R)
       C.append(S_i)
       R = R.difference(S_i)
       costs.append(cost)
       if sum(costs) >= alpha * len(U):
            break 
    return C, costs

def weighted_sample(population, weights, k):
    return random.sample(WeightedPopulation(population, weights), k)

def combination_sample(population, k):
    return itertools.combinations(population, k)

class WeightedPopulation(Sequence):
    def __init__(self, population, weights):
        assert len(population) == len(weights) > 0
        self.population = population
        self.cumweights = []
        cumsum = 0 # compute cumulative weight
        for w in weights:
            cumsum += w   
            self.cumweights.append(cumsum)  
    def __len__(self):
        return self.cumweights[-1]
    def __getitem__(self, i):
        if not 0 <= i < len(self):
            raise IndexError(i)
        return self.population[bisect.bisect(self.cumweights, i)]

def updateSets(Sdict,Unew):
    Sout = {}
    Unew = set(Unew)
    for k,v in Sdict.items():
        nv = set(v).intersection(Unew)
        if nv:
            Sout[k] = nv
    return Sout

# def subSamplingMsmc(mutDict, mutDictInfo , alpha, maxIter = 10000): 
#     ''' U should not be sorted '''
#     out         = defaultdict(int)
#     k           = int(alpha * len(mutDictInfo['mutGintSmpNum'] ) ) 
#     i           = 1 
#     for tempsample  in combination_sample(mutDictInfo['mutGintSmpNum'],k):  
#     # for i in range(numIter):
#         # tempsample     = random.sample(mutDictInfo['mutGintSmpNum'],k) 
#         tempMutDict    = updateSets(mutDict, tempSample)
#         tempResSets, _ = msmc(tempMutDict.values(), set(tempSample))
#         # print tempmutdict.keys(), set(tempsample)
#         # print tempressets
#         tempResGSet  = []
#         for gene, v in tempMutDict.items() :
#             if v in tempressets:
#                 tempResGSet.append(gene)   
#                 tempResStr = ";".join(tempResGSet)
#                 if out.get(tempResStr,''):
#                     out[tempResStr][1] = out[tempResStr][1] + 1
#                 else:
#                     out[tempResStr] = [len(tempResGSet), 1]
#                 print "iteration " + str(i) + " selected genes " + tempResStr    
#             i = i + 1
#             if i > maxIter:
#                 break 
#     out = sorted(out.items(),key=lambda x:x[1][0]) 
#     return out 

def findSsCutUp(alpha, csCutUp):
    i       = 1 
    sseff   = math.ceil(math.sqrt(csCutUp))
    ss      = 2 * sseff 
    while True:
        if scipy.misc.comb(ss, sseff) <= csCutUp:
            return int(ss/alpha)
        else:
            ss = ss - 1 
            i = i + 1
        
def genCombination(population, alpha, ssCutLow, csCutUp):
    ss      = len(population)
    sseff   = int(math.ceil( ss * alpha))
    combs   = scipy.misc.comb( ss, sseff)
    # ssCutup = findSsCutUp(alpha, csCutUp) 
    if sseff <= ssCutLow:
        combSet = population
    elif scipy.misc.comb(ss, sseff) > csCutUp:
        # randPopu = random.sample(population, ssCutup)
        combSet  = random.sample(combination_sample(population,\
                                sseff), csCutUp)
    else:
        combSet = combination_sample(population, sseff) 
    return combSet

#-----------------find all solution
def findAllSol(population, mutDict, \
               alpha = 0.85, ssCutLow = 5, csCutUp = 100000, gssCut = 1, \
               debug = False):
    '''# tempsample    = weighted_sample(mutdictinfo['mutgintsmpnum'], \
                mutdictinfo['mutgintsmpweight'],k) 
        ## weighted sampling methods
    '''
    resultD  = defaultdict(list)
    i = 0 
    for tempsample in genCombination(population,\
                                     alpha, ssCutLow, csCutUp ):
        i = i + 1
        tempmutdict    = updateSets(mutDict, tempsample)
        tempressets, _ = msmc(tempmutdict.values(), set(tempsample))

        tempgeneset = '' 
        cntg        = 0 
        for gene, v in tempmutdict.items() :
            if v in tempressets:
                cntg = cntg + 1
                tempgeneset = tempgeneset+";"+gene   

        tempgval = resultD.get(tempgeneset, '')
        if cntg < gssCut:
            continue 
        if tempgval:
            resultD['tempgeneset'] = [cntg, tempgval[1] + 1] 
        else:
            resultD['tempgeneset'] = [cntg, 1] 
        if debug:
            print tempmutdict.keys(), set(tempsample)
            print tempressets
            print "iteration " + str(i) + " selected " + str(cntg) + \
                    " genes:\t" + tempgeneset
    return resultD

#----optimizing solution

def  getMinSizeSets(resultD):
    while resultD:
        gset    = result.pop()
        ssMin   = 1000
        valMax  = 0 
        gRes    = []
        for k, v in out.items():
            if v[0] < ssMin :
                ssMin = v[0]
                gRes  = k  
            elif v == ssMin:
                gRes.append(k) 
            else:
                continue
        return gRes 

def getMinSizeOverlapSets(result):
    return ''

def selectSets(resultD, type = 'minsize'):
    if   type == "minisize":
        return getMinSizeSets(result)
    elif type == "minsoverlap":
        return getMinSizeOverlapSets(result) 

#------------------test
def __test__():
    alpha = 0.8 
    U = set(range(1,8))
    S = [set([6]), 
         set([1,2]), 
         set([1]), 
         set([1,2]), 
         set([3,4,5]), 
         set([5,6]), 
         set([7])]
    msmc(S,U) 

def __testSample__():
    seq     = ['A', 'B', 'C', 'DDD', 'E']
    seqfreq = [1, 3, 4, 1 ,2] 
    nSample = 100
    alpha   = 0.8 
    k       = int(0.8 * len(seq)) 
    seqSet  = {'1':set(['A','B']),
               '2':set(['B']),
               '3':set(['DDD','A','C']),
               '4':set(['E','A','B']),
               '5':set(['A']),
               '6':set(['C','E'])}
    # seqNew = weighted_sample(seq, seqfreq, k) 
    # seqNew = random.sample(seq, k)
    # print updateSets(seqSet, seqNew) 
    # print findSsCutUp(0.8, 20)
    # print findSsCutUp(0.8, 10000)

    # print "population\t", seq
    # print "all 4 ele subste\t", list(genCombination(seq, 0.8, 2, 100))
    # print "whole set\t", list(genCombination(seq, 0.6, 4, 100))
    # print "get at most 5 set\t", list(genCombination(seq, 0.6, 2, 5))

    print findAllSol(seq, seqSet, ssCutLow = 2, csCutUp = 100, debug = True) 
    print findAllSol(seq, seqSet, alpha = 0.7, ssCutLow = 4, csCutUp = 100, debug = True) 
    print findAllSol(seq, seqSet, alpha = 0.6, ssCutLow = 2, csCutUp = 7, debug = True) 
__testSample__()

