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
usage = 'python ' + sys.argv[0] + ' -s <binding site file>  -a <annotation file> \
-m <mutations> -o <output> '
example = 'python ' + sys.argv[0] + ' -i <input>  -o <output>'
try:
    opts,args = getopt.getopt(argv,"hs:a:m:o:")
except getopt.GetoptError:
    print usage + "\n" + example 
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print usage + "\n" + example 
        sys.exit()
    elif opt in ("-s"):
        bsfile = arg
    elif opt in ("-a"):
        annofile = arg
    elif opt in ("-m"):
        mutfile = arg
    elif opt in ("-o"):
        output = arg
print('Script path:\t'+ sys.argv[0])
print('Input file:\t' + input)
print('Output file:\t'+ output)


annoChrDict = defaultdict(list)
class Gene3UTR:
    def __init__(self, name, chr, utr3s, utr3e):
        self.name = name
        self.chr = chr
        self.utr3s = int(utr3s)
        self.utr3e = int(utr3e)
    def __str__(self):
        return '(%s, %s, %r, %r)' % \
                (self.name, self.chr, self.utr3s, self.utr3e)
class MutRecord:
    def __init__(self, targ, tarchr, tarps, tarpe, mutg, mutps, mutpe, val):
        self.targ = targ
        self.tarchr = tarchr
        self.tarps = tarps
        self.tarpe = tarpe
        self.mutg = mutg 
        self.mutps = mutps
        self.mutpe = mutpe
        self.val = val
    def __str__(self):
        return '(%s, %s, %s, %s)' % (self.targ, self.tarchr, self.tarps, self.mutg) 
# with(open(annofile)) as f:
#     line = f.readline()
#     while line:
#         tgene, tchr, tps, tpe, tstrand =  line.strip().split("\t")
#         tgObj = Gene3UTR(tgene, tchr, tps, tpe )
#         tgObj.strand = tstrand
#         print tgObj
#         annoChrDict[tchr].append(tgObj) 
#         line = f.readline()

mutDict = defaultdict(list)
with(open(mutfile)) as f:
    line = f.readline()
    while line:
        targ, tarchr, tarps, tarpe, mutg, mutps, mutpe, val = \
                line.strip().split("\t",7)
        mutDict[tarchr].append(\
             MutRecord(targ, tarchr, tarps, tarpe, mutg, mutps, mutpe, val) )
        line = f.readline()

cnt = 0
bsDict = {}
with(open(bsfile)) as f:
    line = f.readline()
    line = f.readline()
    while line:
        ptn = r":|\s|\t|\[|,|\]" 
        _, tg, _, _, trbs, trbe, _ = re.split(ptn, line.strip())
        # print tg, trbs, trbe 
        bsDict[tg] = [trbs,trbe]
        cnt = cnt + 1
        if cnt == 10:
            break
        line = f.readline()
         
