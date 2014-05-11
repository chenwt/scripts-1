import numpy as np
class MutSet():
    def __init__(self, gene, mutSample, zscore):
        self.gene = gene
        self.mutSample = mutSample
        self.mutZscore = dict(zip(mutSample,zscore))
    def __str__(self):
        return "(%s, %s)" % (self.gene, self.mutSample)
    def update(self, unionSample):
        for smp in self.mutSample:
            if not smp in unionSample:
                self.mutSample = self.mutSample.remove(smp)
                self.mutZscore = self.mutZscore.remove(smp)
            else:
                pass


def findWegitedMin(S, R, wtype = 'mean'):
    ''' 
    wtype is to select the method to use weight,
    total: the summation of all mutation zscores;
    mean: the mean of all mutation zscores;
    max: the maximization of mutation zscores for each gene
    '''
    ## get the minimum cost set
    minCost = 99999.0
    minElement = -1
    minSet = ''
    minSname = ''
    for i, s in enumerate(S):
        sname_i = s.gene; ss = s.mutSample; sw = s.mutZscore
        ss_i  = set(R).intersection(set(ss))
        if len(ss_i) == 0:
            continue
        sw_i = [sw[a] for a in ss_i ]
        if wtype == 'total':
            cost = 1/reduce(lambda x, y: x+y , sw_i)
        elif wtype == 'mean':
            cost = len(sw_i)/sum(sw_i)
        elif wtype == 'max':
            cost = 1/ reduce(lambda x, y: max(x,y), sw_i)
        else:
            print "ERRor wtype: use default mean; other option total, max or mean"
            cost = len(sw_i)/sum(sw_i)
        if cost < minCost:
            minCost = cost
            minElement = i
            minSname  = sname_i
            minSet  = ss_i 
    return minSname, minSet, minCost

def wgsc(S, U, alpha = 0.8, tol = 0.001, wtype = 'mean'):
    R = U
    C = []
    G = []
    costs = []
    while len(R) != 0:
       g_i, S_i, cost = findWegitedMin(S, R, wtype = wtype)
       C.append(list(S_i))
       G.append(g_i)
       R = list(set(R).difference(set(S_i)))
       costs.append(cost)
       if len(R) <= int((1 - alpha) * len(U)) :
            break
    return G, C, costs

def __test__():
    from collections import defaultdict
    seqSet  = {'G1':['S2','S4','S6'],
                   'G2':['S1','S3'],
                  'G3':['S1'],
                   'G4':['S1'],
                   'G5':['S5'],
                    'G6':['S3']}
    seq     = ['S1', 'S2', 'S3', 'S4', 'S5','S6']
    weightSet = {'G1':[1.0,0.5,1.5],
                    'G2':[2.0, 2.5],
                    'G3':[2.3],
                   'G4':[1.2],
                   'G5':[2.5],
                   'G6':[3.0]}
    setObjL = []
    for sk, ss in seqSet.items():
        setObjL.append(MutSet(sk,ss,weightSet[sk]))
    geneL, smpL, costL = wgsc(setObjL, seq, wtype = "mean")
    geneL, smpL, costL = wgsc(setObjL, seq, wtype = "total")
    geneL, smpL, costL =  wgsc(setObjL, seq, wtype = "max")


