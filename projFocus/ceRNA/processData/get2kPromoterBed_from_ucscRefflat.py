#!/usr/bin/python
#J.HE
'''
input: ucsc downlaoded refflat fasta file
output: bed file for each genes' 3 primer UTR 
purpose : data process
statue: complete
'''

import os, sys, re
from Bio import SeqIO

input , output, kcut = sys.argv[1:4]

kcut = int(kcut)

fhout = open(output, 'w')

for seqObj in SeqIO.parse(input, "fasta"):
    header = seqObj.description.split() 
    gname = header[0].split("_")[-1]
    chr, ps, pe = re.split(r":|-",\
                        header[1].replace("range=", "").replace("chr","")) 
    strand = header[4].replace("strand=", "")

    recOut = "\t".join([chr, str(int(ps) + 2000 - kcut), str(int(pe) + kcut), strand + "|" + gname ])
    fhout.write(recOut + "\n")

fhout.close()
print "[END]"

'''
shell scripts to get unique record
uniq refseq_human_hg19_allGene_promoter1kb_upstrem.bed | awk
'{a[$1"\t"$2"\t"$3]=a[$1"\t"$2"\t"$3]";"$4;}END{for (i in a) {print i"\t"a[i]}}'
|sort -k 1 -k 2n -k 3n |uniq >
refseq_human_hg19_allGene_promoter1kb_upstrem.sorted.uniq.bed
'''
