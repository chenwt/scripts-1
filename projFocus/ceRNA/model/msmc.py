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

def findSsCutUp(alpha, csCutUp):
    i       = 0 
    ss      = int(math.ceil(math.sqrt(csCutUp))) + 1
    sseff   = int(ss * alpha)
    while True:
        # print "combs", scipy.misc.comb(ss, sseff), "size",ss, "effsize", sseff, csCutUp
        if scipy.misc.comb(ss, sseff) <= csCutUp:
            return ss
        else:
            ss      = ss - 1 
            sseff   = int( alpha * ss )
            i = i + 1
        
def genCombination(population, alpha, sseffCutLow, csCutUp):
    ss      = len(population)
    sseff   = int(math.ceil( ss * alpha))
    combs   = scipy.misc.comb(ss, sseff)
    combt   = '2'
    if sseff <= sseffCutLow:
        combSet = [population]
        combt = '1'
    elif combs  > csCutUp and sseff > sseffCutLow:
        combSet = []
        for i in range(csCutUp):
            combSet.append(random.sample(population, sseff))
        combt = '3'
    else:
        combSet = combination_sample(population, sseff) 
    print 'Combination type:\t', combt
    return combSet

#-----------------find all solution
def findAllSol(population, mutDict, \
               alpha = 0.85, ssCutLow = 5, csCutUp = 100000, gssCut = 1, \
               log = False,  debug = False):
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
            if set(v) in tempressets:
                cntg = cntg + 1
                tempgeneset = tempgeneset + gene +";"  

        tempgval = resultD.get(tempgeneset, '')
        if cntg < gssCut:
            continue 
        if tempgval:
            resultD[tempgeneset] = [cntg, tempgval[1] + 1] 
        else:
            resultD[tempgeneset] = [cntg, 1] 
        if debug:
            print "iteration " + str(i) +  "\tpopulation size", len(tempsample), "\tselected " + str(cntg) + \
                    " genes:", tempgeneset
    if log:
        log = {'totalCombCnt': i, 'effss': len(tempsample), "minCost":cntg} 
        return resultD, log
    else:
        return resultD

#----optimizing solution
def integrateGeneSets(gRes, minSize):
    gfinal = {}
    for k in gRes.keys():
        for g in k.split(";"):
            if g and gfinal.get(g,''):
                gfinal[g] = gfinal[g] + 1
            elif g:
                gfinal[g] = 1
            else:
                pass
    gSelect = sorted(gfinal, key = gfinal.get, \
                         reverse = True)[:minSize + 1] 
    return gSelect

def  getMinSizeSets(resultD):
    minSize   = 1000
    gRes    = {} 
    for k, v in resultD.items():
        if v[0] < minSize :
            minSize = v[0]
            gRes    = {} 
            gRes[k] = v  
        elif v[0] == minSize:
            gRes[k] = v  
        else:
            pass
    return integrateGeneSets(gRes, minSize) 

def getMinSizeOverlapSets(result):
    return ''

def getMaxFreq(resultD):
    freqMax = -1 
    minSize = 1000
    gRes    = {} 
    for k, v in resultD.items():
        minSize = min(minSize, v[0])
        if  v[1] > freqMax:
            freqMax = v[1]
            gRes    = {} 
            gRes[k] = v  
        elif v[1] == freqMax:
            gRes[k] = v  
        else:
            pass  
    return integrateGeneSets(gRes, minSize)

def selectSets(resultD, type = 'minsize'):
    if   type == "minsize":
        return getMinSizeSets(resultD)
    elif type == "maxfreq":
        return getMaxFreq(resultD) 
    elif type == "minssOverlap":
        return getMinSizeOverlapSets(resultD) 
    else:
        print "type not match"

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

    # print findSsCutUp(0.8, 40)
    # print findSsCutUp(0.8, 119)

    result  = findAllSol(seq, seqSet,alpha = 0.5,  ssCutLow = 2, csCutUp = 9, \
                         debug = True) 
    # print selectSets(result)
    # result,log  =  findAllSol(seq, seqSet, alpha = 0.6, ssCutLow = 2, csCutUp = \
                              # 100, debug = True, log = True) 
    # print result,log, selectSets(result)
    # result  = findAllSol(seq, seqSet, alpha = 0.6, ssCutLow = 2, csCutUp = 7, debug = False) 
    # print selectSets(result)
    # print selectSets(result, type = "maxfreq")

# __testSample__()

