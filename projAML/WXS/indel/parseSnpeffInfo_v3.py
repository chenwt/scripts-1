#!/usr/bin/env python

"""
  Purpose: Parses INFO field of VCF file to easily human readable tabular format
  Input: Text file in the VCF format
  Output: Tabular format with "parsed.vcf" extension

  Arguments:
  filename - name of table file
  [Y/N] - include samples. Default "Y"

  Note:
  If you wish to change the fields extracted, modify this file
  in the KEYS (fields names), SEP1 and SEP2 to MODIFY section

  To run:
  python parse_info_vcf.py <filename> <Y/N>

"""
import sys
import os.path
import os
import getopt
import gzip


################################################################################
#   July 3, 2011
#   Authors: Vlad Makarov
#   Language: Python
#   OS: UNIX/Linux, MAC OSX
#   Copyright (c) 2012, The Mount Sinai School of Medicine

#   Available under BSD  licence

#   Redistribution and use in source and binary forms, with or without modification,
#   are permitted provided that the following conditions are met:
#
#   Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the following disclaimer.
#
#   Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation and/or
#   other materials provided with the distribution.
#
#   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#   ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#   WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
#   IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
#   INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
#   BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
#   OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
#   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
#   EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
################################################################################


##########################################################################################################

#!/usr/bin/env python


EFFECT='SNPEFF_EFFECT'
#SNPEFF_KEY=Effect_Impact', 'Functional_Class', 'Codon_Change', 'Amino_Acid_change', 'Amino_Acid_length', 'Gene_Name', 'Gene_BioType', 'Coding', 'Transcript', 'Exon', 'GenotypeNum' ]
SNPEFF_KEY=['SNPEFF_IMPACT', 'SNPEFF_FUNCTIONAL_CLASS', 'SNPEFF_CODON_CHANGE', 'SNPEFF_AMINO_ACID_CHANGE', 'SNPEFF_AMINO_ACID_LENGTH', 'SNPEFF_GENE_NAME', 'SNPEFF_GENE_BIOTYPE', 'SNPEFF_CODING', 'SNPEFF_TRANSCRIPT_ID', 'SNPEFF_EXON_ID', 'SNPEFF_GENOTYPENUM' ]


Impacts=['High', 'Moderate', 'Low', 'Modifier' ]
High=['SPLICE_SITE_ACCEPTOR','SPLICE_SITE_DONOR','START_LOST','EXON_DELETED','FRAME_SHIFT','STOP_GAINED','STOP_LOST','RARE_AMINO_ACID']
Moderate=['NON_SYNONYMOUS_CODING','CODON_CHANGE','CODON_INSERTION','CODON_CHANGE_PLUS_CODON_INSERTION','CODON_DELETION','CODON_CHANGE_PLUS_CODON_DELETION','UTR_5_DELETED','UTR_3_DELETED']
Low=['SYNONYMOUS_START','NON_SYNONYMOUS_START','START_GAINED','SYNONYMOUS_CODING','SYNONYMOUS_STOP']
Modifier=['UTR_5_PRIME','UTR_3_PRIME','REGULATION','UPSTREAM','DOWNSTREAM','GENE','TRANSCRIPT','EXON','INTRON_CONSERVED','INTRON','INTRAGENIC','INTERGENIC','INTERGENIC_CONSERVED','NONE','CHROMOSOME','CUSTOM','CDS']


def snnEffParcer(s, leader='(', trailer=')', splitter='|'):
    pk=[]
    if s.find(leader) < 0 or s.find(trailer) < 0:
        return pk
    else:
        start_of_leader = s.index(leader, 0)
        e=s[0:start_of_leader]
        pair= EFFECT + '=' + e
        pk.append(pair)
        #print pair

        end_of_leader = s.index(leader) + len(leader)
        start_of_trailer = s.index(trailer, end_of_leader)
        param = s[end_of_leader:start_of_trailer].split(splitter)

        i=0
        for k in SNPEFF_KEY:
            p=param[i].strip()
            if len(p)>0:
                pair= k + '=' + p
                pk.append(pair)
                #print pair
            i=i+1
        return pk


def getTheMostHigh(lst, leader='(', trailer=')'):
    count=0
    genes_dict = {}
    ranks=[]
    for s in lst:
        start_of_leader = s.index(leader, 0)
        e=s[0:start_of_leader]
        #if e in High:
        #    genes_dict[count]=4
        #elif e in Moderate:
        #    genes_dict[count]=3
        #elif e in Low:
        #    genes_dict[count]=2
        #else:
        #    genes_dict[count]=1
        #count = count+1
        if e in High:
            ranks.append(4)
        elif e in Moderate:
            ranks.append(3)
        elif e in Low:
            ranks.append(2)
        else:
            ranks.append(1)
            
    maxrank=ranks.index(max(ranks))
    return maxrank


###############################  KEYS (fields names) to MODIFY     #######################################

#KEYS=['SNPEFF_IMPACT', 'SNPEFF_EFFECT', 'SNPEFF_FUNCTIONAL_CLASS', 'SNPEFF_CODON_CHANGE',  'SNPEFF_GENE_NAME', 'SNPEFF_GENE_BIOTYPE', 'SNPEFF_CODING', 'SNPEFF_TRANSCRIPT_ID', 'SNPEFF_EXON_ID', 'GMAF']

# SnpEff keys
KEYS=['SNPEFF_GENE_NAME', 'SNPEFF_FUNCTIONAL_CLASS', 'SNPEFF_EFFECT',  'SNPEFF_IMPACT', 'SNPEFF_AMINO_ACID_CHANGE', 'SNPEFF_EXON_ID', 'SNPEFF_CODON_CHANGE']


SEP1=';'
SEP2='='

#############################   NO NEED TO MODIFY BELOW THIS POINT ########################################

def getFh(filename):
  fh = open(filename, "r")
  if filename.endswith('gz'):
    fh = gzip.open(filename, 'rb')
  return fh



""" Converts string to boolean """
def str2bool(v):
  return v.lower() in ["y", "yes", "true", "t", "1"]

""" Helper method to deduplicate the list"""
def dedup(mylist):
    outlist = []
    for element in mylist:
        if element not in outlist:
            outlist.append(element)
    return outlist

""" Parse key-values pairs"""
def parse_field(text, key, sep1, sep2):
    fields = text.strip().split(sep1)
    onefield=set([])
    for f in fields:
        pairs = f.split(sep2)
        if str(pairs[0]) == str(key):
        #if str(pairs[0]).find(str(key)) > -1:
            #op=";".join(dedup(str(pairs[1]).split(',')))
            #onefield.add(op)
            onefield.add(str(pairs[1]))


    if (len(onefield) > 0):
        return ";".join(dedup(onefield))
    else:
        return '.'

def parse_info(filename, include_samples, output_file_base, extout='.txt', comment_char='##'):
    include_samples=str2bool(sys.argv[2])
    outfile = output_file_base+'.' +extout
    fh_out = open(outfile, "w")

    print "Preparing Variant file "
    print "Input file " + filename

    fh=getFh(filename)
    ok=True
    for line in fh:
        line = line.strip()
        if line.startswith(comment_char):
            #fh_out.write(line +'\n')
            print line

        elif line.startswith('CHROM') or line.startswith('#CHROM') :
            fields = line.split('\t')
            samples = []
            sample_names =[]
            if len(fields)<9:
              print ('')
              print "######## ERROR PARSING FILE ##############"
              print "## VCF file format must have at least 9 fields: "
              print "## CHROM POS     ID        REF ALT    QUAL FILTER INFO   FORMAT"
              print "## followed by sample names"
              print "## See http://www.1000genomes.org/wiki/Analysis/vcf4.0 for documentation"
              print "## You may use our parse_info_table.py script to parse INFO in tabular file format "
              print "##########################################"
              ok=False
              break

            if len(fields)>=10:
              samples = fields[9:len(fields)]
              sample_names = '\t'.join(samples).strip()


            header_head='CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER'.strip()
            header_body='\t'.join(KEYS).strip()           
            header_tail='INFO\tFORMAT'.strip()

            if include_samples==True and len(fields) >= 10:
                header_tail=header_tail+ '\t' + sample_names

            header=(header_head+'\t'+header_body+'\t'+header_tail).strip()
            fh_out.write(header.strip()+'\n')

        else:
            fields = line.split('\t')
            if len(fields)>=10:
              samples = fields[9:len(fields)]
              sample_names = '\t'.join(samples).strip()

            chrom = str(fields[0])
            pos = str(fields[1])
            id = str(fields[2])
            ref = str(fields[3])
            alt  = str(fields[4])
            qual = str(fields[5])
            filter = str(fields[6])
            info = str(fields[7])
            format=str(fields[8])
            
            EFF=str(parse_field(info, str('EFF'),SEP1,SEP2))
            EFFS=EFF.split(',')
            index_most_high=getTheMostHigh(EFFS)
            effect_most_high=str(EFFS[index_most_high])

            pk=snnEffParcer(effect_most_high)
            info=info+';' + (';'.join(pk))

            parsed_names=[]
            for k in KEYS:
                parsed_names.append(str(parse_field(info, str(k),SEP1,SEP2)))

            newline_head=chrom+'\t'+pos+'\t'+id+'\t'+ref+'\t'+alt+'\t'+qual+'\t'+filter
            newline_body='\t'.join(parsed_names).strip()
            newline_tail=(info+'\t' +format).strip()
            if include_samples == True and len(fields) >= 10:
                newline_tail=newline_tail+ '\t' + sample_names
            newline = newline_head + '\t' + newline_body + '\t' + newline_tail

            fh_out.write(newline.strip()+'\n')


    fh.close()
    fh_out.close()
    if ok:
        print ('File was saved as ' + outfile)



def main():

  try:
      opts, args = getopt.getopt(sys.argv[1:], "h", ["help"])
  except getopt.error, msg:
      print msg
      print "For help use --help"
      sys.exit(2)
  # process options
  for o, a in opts:
      if o in ("-h", "--help"):
          print __doc__
          sys.exit(0)


  if len (sys.argv) > 2:
    filename=sys.argv[1]
    include_samples=str2bool(sys.argv[2])
    if os.path.exists(filename) and os.path.isfile(filename):
      output_file_base=filename+'.parsed.'
      parse_info(filename=filename, include_samples=include_samples, output_file_base=output_file_base, extout='.var', comment_char='##')
      sys.exit(1)
    else:
      print 'File '  + filename + ' does not exist '
      sys.exit(1)

  else:
    print __doc__
    sys.exit(1)

if __name__ == "__main__":
     sys.exit(main())
