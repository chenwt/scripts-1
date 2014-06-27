#!/usr/bin/python
#J.HE
#Desp: classes and methods for general usages in projceRNA as well as other
# genomic data realted environments 

### reference 
from collections import OrderedDict

class RefGene:
    ''' reference genomic coordinates for each gene, using entrez gene name as
    key '''
    def __init__(self, name, chr, strand, tss, utr3s, utr3e, utr5s, utr5e):
        self.name = name
        self.chr =  chr
        self.strand = strand
        self.tss    = int(tss)
        self.utr3s  = int(utr3s)
        self.utr3e  = int(utr3e)
        self.utr5s  = int(utr5s)
        self.utr5e  = int(utr5e)
    def __hash__(self):
        return hash((self.name, self.tss))
    def __repr__(self):
        return 'RefGene(%s, %r, %r)' % (self.name, self.chr, self.tss)
    def __str__(self):
        return '(%s, %r, %r)' % (self.name, self.chr, self.tss)

class GeneCoord :
    ''' reference genomic coordinates for each gene, using entrez gene name as
    key '''
    def __init__(self, name, chr, tss):
        self.name = name
        self.chr =  chr
        self.tss = int(tss)

    def __hash__(self):
        return hash((self.name, self.tss))
    def __repr__(self):
        return 'Gene(%s, %r, %r)' % (self.name, self.chr, self.tss)
    def __str__(self):
        return '(%s, %r, %r)' % (self.name, self.chr, self.tss)
    def hasCNVInPromoter(self, Cnv, promoterRegionSize = 2000 ):
        if self.chr == Cnv.chr:
            if self.tss + promoterRegionSize < Cnv.ps or \
               self.tss - promoterRegionSize > Cnv.pe : 
                return 0
            else :
                return 1
        else:
            return 0

    def hasCNVIn3utr(self, Cnv):
        if self.chr == Cnv.chr:
            if self.utr3pe < Cnv.ps or self.utr3ps > Cnv.pe : 
                return 0
            else :
                return 1
        else:
            return 0

    def hasCNVIn5utr(self, Cnv):
        if self.chr == Cnv.chr:
            if self.utr5pe < Cnv.ps or self.utr5ps > Cnv.pe : 
                return 0
            else :
                return 1
        else:
            return 0
  
    def hasCNV(self, Cnv, promoterRegionSize = 2000):
        if self.hasCNVInPromoter(self, Cnv, promoterRegionSize = 2000) or \
           self.hasCNVIn3utr(self, Cnv) or\
           self.hasCNVIn5utr(self, Cnv) :
            return 1
        else:
            return 0 
    def hasMutInGene(self, Mutation):
        ''' using binary searching for defaultdict(list) of chrom) '''
        if self.chr == Mutation.chr:
            if self.tss > Mutation.ps and \
               self.tse < Mutation.ps  : 
                return 0
            else :
                return 1
        else:
            return 0


       
    def disp(self):
        print self.name + str(self.chr) + str(self.tss) 

class Mutation :
    def __init__(self, chr, ps):
        self.chr = chr
        self.ps = int(ps)
    def __eq__(self, other):
        return (isinstance(other, self.__class__) \
                and self.__dict__ == other.__dict__)
    def __comp__(self, mut2):
        if self.chr == mut2.chr and self.ps == mut2.ps:
            return 1
        else:
            return 0
    def __repr__(self):
        return 'Mutation(%s,%r)' % (self.chr, self.ps)
    def __str__(self):
        return '(%s,%r)' % (self.chr, self.ps)
    def disp(self):
        print str(self.chr) + "\t" + str(self.ps) 

class CNV :
    def __init__(self, chr, ps, pe, val):
        self.chr = chr
        self.ps = int(ps)
        self.pe = int(pe)
        self.val = float(val)

class Cernet:
    def __init__(self, x):
        self._r = x 
    def __append__(self, x):
        self._r.append(x)
    def getReg(self, gene):
        reg = [gene]
        for r  in self._r:
            if r[0] == gene or r[1] == gene:
               reg.extend(r) 
            else:
                continue
        return list(OrderedDict.fromkeys(reg))          

def chr2Num(chrom):
    if chrom in ['x', 'X'] :
        chrom = 23
    elif chrom in ['y', 'Y']:
        chrom = 24
    elif not chrom in map(str, range(23)):
        chrom = 25
    return int(chrom)

def formatBarcode(barcode):
    ''' TCGA-A1-A0SF-01A-11 to E2.A1L7.01A'''
    return barcode.replace("-",".")[5:16]


# def bSearch(x,list):
#    '''sorted list binary search'''
#     med = list[len(list)/2]
#     l = len(med) 
#     max = list[-1] 
#     if x >= list[med]:
        
