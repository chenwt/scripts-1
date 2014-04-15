import random
import bisect
import collections
from   flattenSeq import flatten
import collections
from   collections import Counter, Sequence


def findMin(S, R):
    S.sort(key=lambda x:len(x),reverse=True)
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
                # print S[minElement], w[minElement]
        except:
            # Division by zero, ignore
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
    # print "Cover: ", C
    # print "Total Cost: ", sum(costs), costs
    return C, costs

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


def weighted_sample(population, weights, k):
    return random.sample(WeightedPopulation(population, weights), k)

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

def __testSample__():
    seq     = ['A', 'B', 'C', 'DDD', 'E']
    seqfreq = [1, 3, 4, 1 ,2] 
    nSample = 100
    alpha   = 0.8 
    k       = int(0.8 * len(seq)) 
    seqSet  = {'1':['A','B'],
               '2':['B'],
               '3':['DDD','A','C'],
               '4':['E','A','B'],
               '5':['A'],
               '6':['C','E'] }
    # seqNew = weighted_sample(seq, seqfreq, k) 
    seqNew = random.sample(seq, k)
    print seqNew 
    print updateSets(seqSet, seqNew) 

# __testSample__()


