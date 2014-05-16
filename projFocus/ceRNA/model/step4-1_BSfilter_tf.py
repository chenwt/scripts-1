#!/usr/bin/python
#J.HE
# input : 1> the predicted binding site genomic coordinates, 
#        2> somatic mutation site genomic coordinate
# output: 1> the mutations occurs in binding site, and the target gene, as well as the modulator(TF/miRNA)
#Desp.: given predicted binding site, find the mutations that are in the binding
#sites for each gene

import sys,getopt
import re
from collections import defaultdict
from generalUtils import *

# argv = sys.argv[1:]
# input = ''
# output = ''
# usage = 'python ' + sys.argv[0] + '\
#         -s <binding site file> \
#         -m <mutations> \
#         -o <output> '
# example = 'python ' + sys.argv[0] + ' -i <input>  -o <output>'

# try:
#     opts,args = getopt.getopt(argv,"hs:m:o:")
# except getopt.GetoptError:
#     print usage + "\n" + example
#     sys.exit(2)
# for opt, arg in opts:
#     if opt == '-h':
#         print usage + "\n" + example
#         sys.exit()
#     elif opt in ("-s"):
#         tfBSfile = arg
#     elif opt in ("-m"):
#         mutfile = arg
#     elif opt in ("-o"):
#         outputfile = arg

# print('Input file:\t' + tfBSfile)
# print('Input file:\t' + mutfile)
# print('Output file:\t'+ outputfile)

##test file
mutfile="/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/03102014/tcgal2som/genome.wustl.edu__Illumina_All.maf"
tfBSfile="/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/tfBS/TFtarg_binding_site.hg19"
outputfile="/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/fucFilt/mut_in_TFBindSite_05142014.hg19"
##

class GeneRegion:
    def __init__(self, gene, chr, ps, pe, strand):
        self.gene = gene
        self.chr = chr
        self.ps = int(ps)
        self.pe = int(pe)
        self.strand = strand
    def __eq__(self,other):
        return self.gene == other.gene and \
                self.chr == other.chr and \
                self.ps == other.ps and  \
                self.pe == other.pe
                
    def __ne__(self, other):
        return self.gene == other.gene or \
                self.chr == other.chr or \
                self.ps == other.ps or \
                self.pe == other.pe
     
    def __str__(self):
        return '%s\t%s:%s-%s' % (self.gene, self.chr, self.ps, self.pe )
    
    def __expr__(self):
        return '%s\t%s:%s-%s' % (self.gene, self.chr, self.ps, self.pe )
    
    def __cmp__(self, grObj):
        if self.chr == grObj.chr:
            if self.pe < grObj.ps:
                return -1
            elif self.ps > grObj.pe:
                return 1
            else :
                return 0
        else:
            return -2
    
    def __intersect__(self, gRObj):
        if self.chr == gRObj.chr:
            if self.pe < gRObj.ps or self.ps > gRObj.pe:
                return 0
            else:
                return 1
        else:
            return 0
    def info(self):
        return(self.gene + "\t" + self.chr + ":" + "-".join(map(str),[self.ps, self.pe]))

def __test_GRObj():
    testGG1 = GeneRegion('ZYX', 1, 23, 23 + 1, '+')
    testGG2 = GeneRegion('X', 1, 23, 30, '+')
    testGG3 = GeneRegion('Z', 1, 10, 22,'-')
    testGG4 = GeneRegion('X', 1, 23, 30, '+')
    
    print testGG1 == testGG2
    print testGG2 == testGG4
    print testGG1.__intersect__(testGG2)
    print testGG1.__intersect__(testGG3)

###load mutation info
mutDict = defaultdict(list)
cnt = 0
with(open(mutfile)) as f:
    head = f.readline()
    line = f.readline()
    while line:
        g_crt, _, _,  _, chr_crt, ps_crt, pe_crt, strand_crt = \
                line.strip().split("\t",7)
        crtMut = GeneRegion(g_crt, chr_crt, ps_crt, pe_crt, strand_crt)
        mutDict[g_crt].append(crtMut)
        cnt = cnt + 1
        line = f.readline()
print "mutations loaded : %s" % cnt


def writeLine(fileH, grObj1, grObj2, mod):
    fileH.write( grObj1.gene + "\t" + grObj1.chr + \
                ":" + "-".join(map(str,[grObj1.ps, grObj1.pe])) + \
                "\t" + grObj2.gene + "\t" + grObj2.chr + \
                ":" + "-".join(map(str,[grObj2.ps, grObj2.pe])) + \
                "\t" + mod + "\n" ) 

bsDict = defaultdict(list)

outputH = open(outputfile, 'wt')

cnt = 0
mutOrigDict = mutDict
cntFind = 0 

with(open(tfBSfile)) as f:
    modType = f.readline().strip().split("\t")[0]
    outHeader = ["mutGene", "mutChr:mutPs-mutPe","targetGene","bsChr:bsPs-bsPe", "modulator" + modType]
    outputH.write("\t".join(outHeader) + "\n" )
    line = f.readline()
    while line:
        mod_crt, g_crt, coor_crt, _ = line.strip().split("\t", 3)
        chr_crt, ps_crt, pe_crt = re.split(":|-", coor_crt.replace('chr',''))
        bsObj_crt = GeneRegion(g_crt, chr_crt, ps_crt, pe_crt,'+')
        bsLen = abs(int(pe_crt) - int(pe_crt))
        for crtmut in mutDict[g_crt]:
            if crtmut.__intersect__(bsObj_crt) :
                writeLine(outputH, crtmut, bsObj_crt, mod_crt)
                mutDict[g_crt] = filter(lambda x:x.__eq__(crtmut), mutDict[g_crt])
                cntFind  = cntFind + 1
        cnt  = cnt + 1
        line = f.readline()

print "mutations find in tf binding area:\t %s" % cntFind
outputH.close()
print "#[END]"
