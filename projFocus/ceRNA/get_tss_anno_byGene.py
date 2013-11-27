#!/usr/bin/env python
#input:<file: entrez_tss.bed>
#output:<file: entre_tss_by_gene.bed>
#desp. concantenate tss for each gene

import os
import sys

argv = sys.argv[1:]
in_file = argv[0]
out_file = in_file + ".signleTSS"
#in_file = "Transcription_Start_Sites_July_2010_hg19.gff_EntrezID.bed"
#out_file = "Transcription_Start_Sites_July_2010_hg19.gff_EntrezID_by_Gene.bed"
#in_file = "test.txt"
#out_file ="test_out.txt"

outf = open(out_file,"w")
d = {}
with open(in_file) as inf1:
  for i, line in enumerate(inf1):
    if i == 0:
      outf.write("Gene" + "\t" + "chrom" +"\t" + "start" +"\t"+ "strand" + "\n")
    else: 
      [key, chrom, pos,strand] = line.strip().split("\t")
      if key in d.keys():
        [chrom_prev,pos_prev,strand_prev] = d[key]
        if chrom_prev == chrom :
          if (strand == "+" and int(pos) < int(pos_prev)) or ( strand == "-" and int(pos) > int(pos_prev)): 
            d[key] = [chrom, pos, strand]
          else:
            continue
        else:
          continue
      else :
        d[key] = [chrom, pos, strand]
  
for key in sorted(d.keys()):
  outf.write(key + "\t" + "\t".join(d[key]) + "\n")




