class MutRecord:
    '''
    Mutation Records Objects, store target gene, chr, mutations gene, mutation
    starting site ,mutation end site and mutatinon spectrum of all avaialbe
    samples 
    more details
    '''
    def __init__(self, targ, tarchr,  mutg, mutps, mutpe, val):
        self.targ = targ
        self.chr = tarchr
        # self.tarps = tarps
        # self.tarpe = tarpe
        self.mutg = mutg 
        self.ps = int(mutps)
        self.pe = int(mutpe)
        self.val = val

    def __eq__(self,other):
        return self.targ == other.targ and self.chr == other.chr \
                and self.ps == other.ps and self.pe == other.pe
        
    def updateVallist(self, mut2):
        ##make sure valarray is np array 
        return self.valarray + mut2.valarray
    def __str__(self):
        return '(%s, %s, %s, %s)' % (self.mutg, self.chr, self.ps, self.pe ) 
    def inRegion(self, regObj):
        if self.chr == regObj.chr:
            if self.pe < regObj.ps or self.ps > regObj.pe:
                return 0
            else:
                return 1 
        else:
            return 0 
    def allInfo(self):
        return "\t".join(map(str, [self.targ, self.mutg, self.chr, self.ps, \
                                   self.pe] )) + "\t" + self.val 


