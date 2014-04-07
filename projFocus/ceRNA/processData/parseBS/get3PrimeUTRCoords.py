#!/usr/bin/python
#J.HE
'''Descp.: Given binding site sequence infor, search the genomic coordinaites in
DNA sequanece, get the coordinate of miRNA binding sites for each gene'''

import sys,getopt
import re
from collections import defaultdict
from searchStr import bruteSearch 
argv = sys.argv[1:]
input = ''
output = ''
usage = ""
example = ""
try:
    opts,args = getopt.getopt(argv,"hc:d:o:")
except getopt.GetoptError:
    print usage + "\n" + example 
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print usage + "\n" + example 
        sys.exit()
    elif opt in ("-c"):
        cupidseq = arg
    elif opt in ("-d"):
        dnaseq = arg
    elif opt in ("-o"):
        output = arg
print('Script path:\t'+ sys.argv[0])
print('Input file:\t' + cupidseq)
print('Input file:\t'+ dnaseq)
print('Output file:\t'+ output )

##load all cupidseq 
bsSeqDict = defaultdict(list)
with(open(cupidseq)) as f:
    line = f.readline()
    line = f.readline()
    while line:
        gene, bsseqinfo = line.strip().split("\t")
        for x in bsseqinfo.split(";"):
            bsSeqDict[gene].append(x.lower())
        line = f.readline()
print "binding seq loaded"
# print bsSeqDict.items()[1]

##process DNA seq by gene
def find_all(qstr, allstr):
    start=0
    while True:
        start = allstr.find(qstr, start)
        if start == -1: return
        yield start
        start += len(qstr)

outputH = open(output, 'w')
outputH.write("Symbol\tChr:bsStart-bsEnd\n")
cnt = 0 
with(open(dnaseq)) as f:
    line = f.readline()
    while line:
        cnt = cnt + 1
        if cnt % 1000 == 0 :
            print " %s line processed" % cnt
        if not re.match("^Symbol",line):
            gene, coord, seq  = line.strip().split("\t")
            if bsSeqDict.get(gene, ''):
                outbss = []
                for bsseq in bsSeqDict[gene]:
                    bsstart, querySeq =  bsseq.split(":")
                    bsstart = int(bsstart)
                    for bsindex in find_all(querySeq, seq): 
                        # print querySeq
                        # print seq[bsindex:bsindex+57]
                        chrom, tss, tse = re.split(":|-", coord) 
                        outbss.append(int(tss) + bsindex + 1 - bsstart)
                outRec = gene + "\t" +  chrom + ":" + \
                        ";".join(map(str, list(set(outbss))))
                outputH.write(outRec + "\n" )
                del bsSeqDict[gene]
        line = f.readline()
outputH.close()
