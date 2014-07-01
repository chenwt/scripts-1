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

input , output = sys.argv[1:3]

fhout = open(output, 'w')

for seqObj in SeqIO.parse(input, "fasta"):
    header = seqObj.description.split() 
    gname = header[0].split("_")[-1]
    chr, ps, pe = re.split(r":|-",\
                        header[1].replace("range=", "").replace("chr","")) 
    strand = header[4].replace("strand=", "")
    recOut = "\t".join([chr, ps, pe, strand + "|" + gname ])
    fhout.write(recOut + "\n")
    # if int(ps) + 30 < int(pe):
        # recOut = "\t".join([chr, str(int(ps) + 30), pe, strand + "|" + gname ])
        # fhout.write(recOut + "\n")

fhout.close()
print "[END]"

