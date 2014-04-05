#!/usr/bin/python
#J.HE
'''Descp.: Given binding site sequence infor, search the genomic coordinaites in
DNA sequanece, get the coordinate of miRNA binding sites for each gene'''

import sys,getopt
from collections import defaultdict
from searchStr import bruteSearch 
argv = sys.argv[1:]
input = ''
output = ''
try:
    opts,args = getopt.getopt(argv,"hm:d:o:")
except getopt.GetoptError:
    print usage + "\n" + example 
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print usage + "\n" + example 
        sys.exit()
    elif opt in ("-m"):
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
        bsSeqDict[gene].append(bsseqinfo.split(";")) 
        line = f.readline()
##process DNA seq by gene

with(open(dnaseq)) as f:
    line = f.readline()
    while line:
        gene, coord, seq  = line.strip().split"\t")
        seq =  seq.upper()
        if bsSeqDict.get(gene, ''):
            for bsseq in bsSeqDict[gene]:
                bsstart, querySeq =  bsseq.split(":")
                bsindex = bruteSearch( querySeq, seq) 
                if bscoord:
                    chrom, tss, tse = re.split(":|-", coord) 
                    tss = int(tss)
                    print gene + "\t" +  chrom + ":" + str(tss + bsindex + 1)\
                            + "-" + str(tss+bsindex + 1 + 7) 

        line = f.readline()
