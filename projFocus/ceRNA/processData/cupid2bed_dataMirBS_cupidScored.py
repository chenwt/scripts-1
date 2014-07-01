#!/usr/bin/python
#J.HE
'''
    input: 1. cupid scored prediction in $fdata/miBS, 
       2. 3'UTR RNA sequence(57 bp) , 
       3. CLASH gene 3UTR DNA fasta file, 
    output: bed like file of binding site coordinates, including UTR lenght and score
    purpose: data process, under work
    last modified: Jun 28, 2014
'''

import os, sys, re
from collections import defaultdict

argv = sys.argv[1:]

input = argv[0]
seqRna  = argv[1]
seqDna  = argv[2]

output = argv[-1]


## step 1. get gene binding sequence, file 1, file 2, 

##--store all 3prime utr sequence into memory, require memory
rnaseqDict = defaultdict(list)
with(open(seqRna)) as f:
    line = f.readline()
    if re.findall(r"(?i)symbol", line):
        line = f.readline()
    while line:
        gname, _, _, rnaseq = line.strip().replace("Entrez:",\
        "").replace("HGNC:","").split("\t")
        rnaseqDict[gname].append(rnaseq) 
        line = f.readline()


def findBSseq( DNAseq, bss_bse):
    ''' function to extract 57 bp sequnce from utr seq 
        return the binding site in the flanking seq, and the 57bp sequence
    '''
    bss, bse = map(int, bss_bse.split("-"))
    ldna = len(DNAseq) 
    if ldna < 57 or ldna <= bss + 7 :
        # print "tot len fail"
        return "", "" 
    elif bss <= 25:
        # print "bss lt 25"
        outseq = DNAseq[:57]
        return str(bss) + "-" + str(bss + 7 -1), outseq
    elif ldna - bse < 25 and ldna > bss + 7:
        # print "bse lt 25 to the end"
        outseq = DNAseq[-57:] 
        return  str(57 - ( ldna - bse) - 6) + "-" + str(57 - (ldna - bse) ), outseq 
    else:
        # print "in the mid"
        outseq =  DNAseq[bss-1-25: bse + 25] 
        return str(25) + "-" + str(25 + 7) , outseq 

###----load  binding site infomration, extract binding sequence, write temp file 
tempfile = input + ".temp"
fhtemp1 = open(tempfile, 'w') 
cnt = 0
with(open(input)) as f:
    line = f.readline()
    if re.findall(r"(?i)gene", line) or re.findall(r"(?i)symbol",line):
        line = f.readline()
    tempOutBlock = [] 
    while line:
        score, gname, refseqid, mirid, bs = \
                line.strip().replace("HGNC:","").replace("Entrez:","").split()

        if not rnaseqDict.get(gname, 0) :
            line = f.readline()
            continue 
        bs = bs.replace("[","").replace("]","").replace(",", "-")
        for v in rnaseqDict[gname]: 
            bsinseq, bsseq = findBSseq(v, bs) 
            if bsinseq and max(map(int, bsinseq.split("-"))) > 57:
                print gname, refseqid, score, mirid, bs, bsinseq 
                exit()
                continue 
            tempOutBlock.append("\t".join([gname, refseqid, score, mirid,\
                                           bs, bsinseq, bsseq]))
            if len(tempOutBlock) == 100 : 
                # print tempOutBlock
                fhtemp1.write( "\n".join(tempOutBlock) + "\n" )
                tempOutBlock = [] 
                cnt = cnt + 100
    
        line = f.readline()
fhtemp1.write( "\n".join(tempOutBlock) + "\n" )
fhtemp1.close()

cnt = cnt + len(tempOutBlock)
print "write temp file record", cnt
del tempOutBlock
del cnt 

###--release rnaseq memory and load bs seq into memory
del rnaseqDict

bsseqDict = defaultdict(list)
bsinfodict = defaultdict(list)
bsPosInseqDict =  defaultdict(list)
with(open(tempfile)) as f:
    line = f.readline()
    while line:
        line =  line.strip().split("\t")

        if len(line) < 7:
            line = f.readline()
            continue
        gname = line[0]; score = line[2]; 
        mirid = line[3]; bsinseq = line[5]
        bsseq = line[6];
    
        bsinfodict[gname].append(score+"\t"+mirid)
        bsPosInseqDict[gname].append(bsinseq)
        bsseqDict[gname].append(bsseq)
        line = f.readline()

print "bs seq temp file loaded"

###---parse fasta file
from Bio import SeqIO 
from searchStr import searchUTR

fhout = open(output, 'w')
# outHead = "\t".join(["chr", "psBS", "peBS", "strand", "geneName","score","mirId"])
# fhout.write(outHead + "\n")
cnt = 0 
for seqRec in SeqIO.parse(seqDna, 'fasta'):
    header = seqRec.description.split()

    gname = header[0].split("_")[2]
    if not bsseqDict.get(gname, 0):
       continue  
        
    chr, ps, pe = re.split("-|:", header[1].replace("range=chr",""))
    ps = int(ps)
    pe = int(pe)
    strand = header[4].replace("strand=","")

    ## get rna seq match pos in dna seq
    for (rnas, rnas_bspos, rnas_info) in \
                zip(bsseqDict[gname], bsPosInseqDict[gname], bsinfodict[gname]):
        for posbs in searchUTR(rnas, str(seqRec.seq)) :
            inseqstart, inseqend  = map(int, rnas_bspos.split("-"))
            outRec = [chr, ps + posbs  + inseqstart, ps + posbs + inseqend, \
                strand, gname, rnas_info ] 
            cnt = cnt + 1
            fhout.write("\t".join(map(str,outRec) ) + "\n")

fhout.close() 
print "output record:" , cnt 
print "[END]"
