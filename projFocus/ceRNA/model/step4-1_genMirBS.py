#!/usr/bin/python
#J.HE
# given cupid predicted microRNA binding sites,
# download 3'UTR sequence,
# generate the predicted genomic coordinate of binding sites

import os, sys, getopt
import re
from collections import defaultdict
from searchStr import bruteSearch
from searchStr import searchUTR

argv = sys.argv[1:]
input = ''
output = ''
usage = 'python ' + sys.argv[0] + ' -i <input>  -o <output>'
example = 'python ' + sys.argv[0] + ' -i <input>  -o <output>'
try:
    opts,args = getopt.getopt(argv,"hc:u:d:o:")
except getopt.GetoptError:
    print usage + "\n" + example
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print usage + "\n" + example
        sys.exit()
    elif opt in ("-c"):
        cupidfile = arg
    elif opt in ("-u"):
        seqfile = arg
    elif opt in ("-d"):
        dnaseqfile = arg
    elif opt in ("-o"):
        outputfile = arg
        
print('Input file--------------------:')
print(cupidfile)
print(seqfile)
print(dnaseqfile)
print('Output file---------')
print(outputfile)


##-test
# seqfile = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/data/cupidPrediction/3PrimeUTR"
# cupidfile = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/data/cupidPrediction/CupidPred"
# dnaseqfile = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/data/cupidPrediction/tempfile/refseq_hg19_refflat_cupid3pGene.fasta_Apr4.tsv_1"
# outputfile = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/data/miBS/mircoRNA_BindSite_cupidPredict.hg19.test"

##load 3'UTR seqinfor
seq3pUTR = defaultdict(list)
cnt = 0
with(open(seqfile)) as f:
    line = f.readline()
    while line:
        if not re.match("^Symbol", line):
            _, genename, refseqid, _, seq  = re.split(r":|\s+",line.strip())
            seq3pUTR[genename].append(seq)
            cnt = cnt+1
        line = f.readline()
print "%s sequence for %s genes" % (cnt, len(seq3pUTR.keys()) )


##get cupid prediction 57 bp seq
miReBSD = defaultdict(list)

cnt = 0
with(open(cupidfile)) as f:
    head = f.readline()
    line = f.readline()
    while line:
        _, genename, refseqid, mir_crt, ps_crt, pe_crt = \
                re.split(r":|\s+|,", line.strip().replace("[","").replace("]",""))
        ps_crt = int(ps_crt)
        pe_crt = int(pe_crt)
        for gSeq in  seq3pUTR[genename]:
            ps = max(ps_crt  - 1 - 25, 0)
            pe = min(pe_crt + 25, len(gSeq) )
            str(ps_crt - ps) + ":" + gSeq[ps:pe] 
#             print genename, refseqid, mir_crt, ps_crt, pe_crt, ps, pe, gSeq[ps:pe]
            miReBSD[genename].append([mir_crt, ps_crt - ps, pe_crt - ps , gSeq[ps:pe]])
            cnt = cnt + 1
        line = f.readline()


##process DNA seq by gene

def writeline(fileH, gene, mir, chr, ps, pe):
    outRecord = mir + "\t" + gene + "\t" + \
                str(chr) + ":" + str(ps) + "-" + str(pe)
    fileH.write(outRecord + "\n")

outH = open(outputfile, 'wt')
writeline(outH, "targetGene", "microRNA", "chr", "bindStart", "bindEnd")

cnt = 0
cntout = 0
with(open(dnaseqfile)) as f:
    line = f.readline()
    while line:
        cnt = cnt + 1
        if cnt % 100 == 0 :
            print " %s line processed" % cnt
        gene, coord, seq  = line.strip().split("\t")
        chrom, tss, tse = re.split(":|-", coord)
        tss = int(tss)
        miBSseq_crt = miReBSD.get(gene, '')
        if miBSseq_crt:
            for mir, bps, bpe, querySeq in miBSseq_crt:
                for bsindex in searchUTR(querySeq, seq):
                    coord_bps = tss + bsindex  + bps 
                    coord_bpe = tss + bsindex  + bpe
                    writeline(outH, gene, mir, chrom, coord_bps, coord_bpe)
                    cntout = cntout + 1
            del miReBSD[gene]
        line = f.readline()
      
outH.close()

print "mir-target pair" , cntout
print "#[END]"


