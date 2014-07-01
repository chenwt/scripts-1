#!/usr/bin/python
#J.HE
'''
input: 1 CLASH_6102 in $fdata/miBS, 
       2. 3'UTR RNA sequence(57 bp), 
       3. CLASH gene 3UTR DNA fasta file, 
output: bed like file of binding site coordinates
purpose: data process, under work
last modified: Jun 27, 2014
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

###----load  binding site infomration, extract binding sequence, write temp file 
fhtemp1 = open("temp1", 'w') 
cnt = 0
with(open(input)) as f:
    line = f.readline()
    if re.findall(r"(?i)gene", line):
        line = f.readline()
    tempOutBlock = [] 
    while line:
        
        gname, refseqid, len3putr, bs = \
                line.strip().replace("HGNC:","").replace("Entrez:","").split("\t")
        if not rnaseqDict.get(gname, 0) :
            line = f.readline()
            continue 

        for v in rnaseqDict[gname]: 
            bss, bse = map(int, bs.split("-"))
            tempOutBlock.append("\t".join([gname, refseqid, len3putr,\
                                           bs, v[bss-1:bse]]))
        if len(tempOutBlock) == 20 : 
            fhtemp1.write( "\n" + tempOutBlock + "\n" )
            tempOutBlock = [] 
            cnt = cnt + 20
    
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

with(open("temp1")) as f:
    line = f.readline()
    while line:
        line =  line.strip().split("\t")

        if len(line) < 5:
            line = f.readline()
            continue
            
        gname = line[0]; bsseq = line[4]; bs = line[3]

        bsinfodict[gname].append(bs)
        bsseqDict[gname].append(bsseq)
        line = f.readline()

print "bs seq temp file loaded"

###---parse fasta file
from Bio import SeqIO 
from searchStr import searchUTR

fhout = open(output, 'w')
outHead = "\t".join(["chr", "psBS", "peBS", "strand"])
fhout.write(outHead + "\n")
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
    for rnas in bsseqDict[gname]:
        for posbs in searchUTR(rnas, str(seqRec.seq)) :
            outRec = [chr, ps + posbs, ps + posbs + len(rnas) - 1, \
                strand, gname] 
            fhout.write("\t".join(map(str,outRec) ) + "\n")
fhout.close() 
print "[END]"
