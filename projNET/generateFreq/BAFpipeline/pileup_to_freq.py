#! /usr/bin/env python
# pileup_to_freq.py is called by pileup_to_freq.sh.  Execute ./pileup_to_freq.sh for documentation for this file.   

import geneutils
from optparse import OptionParser

def MakeFreq(pileup_file, output_file):

  PILEUP_ob = geneutils.PILEUPFileOb(memory_flag = False, fn = pileup_file)

  SNP_list = []
  for v in PILEUP_ob:
    SNP_list.append([v.getPos_t()[0], v.getPos_t()[1], v.getReadCount(),\
                    v.getAltNucleotideCount()])
  summary = ('@this file ', output_file,\
             ' was created by make_freq.py on ',\
             geneutils.getTimeStamp(), ' from the pileup file ',\
             pileup_file, '\n') 
  with open(output_file, 'w') as fh:
    fh.write("".join(summary))
    fh.write("@\n@\n@\n@\nChr\tPos\tTotalreads\tAltreads\n")
  geneutils.appendListToFile(SNP_list, output_file)

if __name__ == '__main__':
  parser = OptionParser()
  (options, args) = parser.parse_args()
  pileup_file = args[0]
  output_file = args[1]  
  MakeFreq(pileup_file, output_file)  

