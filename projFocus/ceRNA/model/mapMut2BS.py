#!/usr/bin/python
#J.HE
#Desp.: given predicted binding site, find the mutations that are in the binding
# sites for each gene


import sys,getopt
import re
# import textwrap
# print textwrap.fill(usage, 70)
from collections import defaultdict
from generalUtils import GeneCoord

argv = sys.argv[1:]
input = ''
output = ''
usage = 'python ' + sys.argv[0] + ' -s <binding site file> \
-m <mutations> -o <output> '
example = 'python ' + sys.argv[0] + ' -i <input>  -o <output>'
try:
    opts,args = getopt.getopt(argv,"hs:m:o:")
except getopt.GetoptError:
    print usage + "\n" + example 
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print usage + "\n" + example 
        sys.exit()
    elif opt in ("-s"):
        bsfile = arg
    elif opt in ("-m"):
        mutfile = arg
    elif opt in ("-o"):
        output = arg
print('Script path:\t'+ sys.argv[0])
print('Input file:\t' + bsfile)
print('Input file:\t' + mutfile)
print('Output file:\t'+ output)


class MutRecord:
    # def __init__(self, targ, tarchr, tarps, tarpe, mutg, mutps, mutpe, val):
    def __init__(self,   targ, tarchr,  mutg, mutps, mutpe, val):
        self.targ = targ
        self.chr = tarchr
        # self.tarps = tarps
        # self.tarpe = tarpe
        self.mutg = mutg 
        self.ps = int(mutps)
        self.pe = int(mutpe)
        self.val = val
    def __eq__(self,other):
        return self.mutg == other.mutg and self.chr == other.chr \
                and self.ps == other.ps and self.pe == other.pe
        
    def __str__(self):
        return '(%s, %s, %s, %s)' % (self.tarchr, self.mutg, self.ps, self.pe ) 
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

class Region:
    def __init__(self, chr, ps, pe):
        self.chr = chr
        self.ps = int(ps)
        self.pe = int(pe)

mutDict = defaultdict(list)
with(open(mutfile)) as f:
    head = f.readline()
    line = f.readline()
    while line:
        # targ, tarchr, tarps, tarpe, mutg, mutps, mutpe, val = \
        targ, tarchr, _,  _, mutg, mutps, mutpe, val = \
                line.strip().split("\t",7)
        crtMut = MutRecord(targ, tarchr,  mutg, mutps, mutpe, val) 
        if not crtMut in mutDict[mutg]: 
            mutDict[mutg].append(crtMut)
        line = f.readline()
print "mutation matrix loaded : %s" % len(mutDict.keys())
bsDict = defaultdict(list)
outputH = open(output, 'wt')
outputH.write("\t".join(\
            ["geneTar", "geneMut", "chrome", "psMut", "peMut"] ) \
            + "\t" +  head.split("\t",7)[7]   )
cnt = 0 
with(open(bsfile)) as f:
    line = f.readline()
    while line:
        if not line.strip().endswith(':') :
            gene, chrom, bsPos = re.split(":|\t", line.strip())
            for crtmut in mutDict[gene]:
                for x in bsPos.split(";"):
                    x = int(x)
                    if crtmut.inRegion(Region(chrom, x, x + 7)) :
                        outputH.write(crtmut.allInfo() + "\n") 
                        cnt = cnt + 1
                        break
        line = f.readline()

print "mutations find in binding area:\t %s" % cnt 
outputH.close()

