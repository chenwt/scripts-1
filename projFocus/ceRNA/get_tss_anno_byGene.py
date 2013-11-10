#!/usr/bin/env python
#input:<file: entrez_tss.bed>
#output:<file: entre_tss_by_gene.bed>
#desp. concantenate tss for each gene

in_file = "Transcription_Start_Sites_July_2010_hg19.gff_EntrezID.bed"
out_file = "Transcription_Start_Sites_July_2010_hg19.gff_EntrezID_by_Gene.bed"
#in_file = "test.txt"
#out_file ="test_out.txt"

outf = open(out_file,"w")
d = {}
with open(in_file) as inf1:
  for i, line in enumerate(inf1):
    [key, val1, val2, val3] = line.strip().split("\t")
    if key in d.keys():
      [v1,v2,v3] = d[key]
      v2 = v2 + ":" + val2
      d[key] = [v1,v2,v3]
      #print key, "\t".join(d[key])
      outf.write(key + "\t" + "\t".join(d[key]) + "\n")
    else :
      d[key] = [val1, val2, val3]
      #print key,"\t".join(d[key])
      outf.write(key + "\t" + "\t".join(d[key]) + "\n")




