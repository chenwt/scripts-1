#!/usr/bin/env python
#input: 1.<file:vcf after gatk,annovar, single calling>
#	2.<file:potential false positive gene names, one by one line>
#	3.<string: output file name> 
#output: 1.<file: maf after filtering> 
#desp: filtering each vcfs and transform it into maf(acceptable by mutSig)
##TODO: 1.do filtering at first step, do maf in another file 


#loading package
import os
import sys,getopt
import re
import collections

#get input
# #-------test_begin-----
# inp = 'input_test.vcf'
# out = 'out_test.tsv'
# outlog = out + ".log"
# #------test_end-------

usage="[UASGE:]" +sys.argv[0] + " -i <input.vcf> -f <falsePostivegene.txt> -o <out.maf>"
err="[ERR:]"
example="[EXAMPLE:]" + sys.argv[0] + " -i -f -o "
argv = sys.argv[1:]
try:
      # opts,args = getopt.getopt(argv,"hi:o:f:")
      opts,args = getopt.getopt(argv,"hi:o:")
except getopt.GetoptError:
      print usage + "\n" + example 

for opt,arg in opts:
      if opt == "-h":
	    print usage + "\n" + example 
      elif opt == "-i":
            inp = arg
      # elif opt == "-f":
            # fpgenefile = arg
      elif opt == "-o":
            out = arg

#-----------------setting parameter
##---genomic level
pltum_cut = 40
clr_cut = 1
dp3_cut = 1
dp4_cut = 1
dp34_cut = 3
dp_cut = 10
mq0_cut =  0.1
# gtnorm_cut = '0/0'
# plnorm_cut = 30
# dpnorm_cut = 1
# mafnorm_cut = 0.05
maftum_cut = 0
##-----variant level
dbsnp = 'DB'
variantType = 'SNP'
funcClass_type = 'synonymousSNV'
snv_type = ''
##----knowledge base level

fpgenefile = "/ifs/data/c2b2/ac_lab/jh3283/projFocus/other/positiveControlGene.txt" 
fpgeneArray = []
with open(fpgenefile) as fpgenef:
  for line in fpgenef.readlines():
    fpgeneArray.append(line.strip())
print "fpgene Num: " + str(len(fpgeneArray))

# mafHeader=[]
# with open(mafHeaderFile) as mfH:
#   for line in mfH.readlines():
#     mafHeader = line.strip().split("\t") 
# print ",".join(mafHeader)

outf = open(out,'w')
###---.maf file,
'''
Hugo_Symbol
Tumor_Sample_Barcode
effect (as one of ):nonsilent /  silent / noncoding
categ (one of ): CpG transitions/  CpG transversions/ C:G transitions/  C:G transversions/ A:T transitions/A:T transversions/null+indel mutations
 Variant_Classification, Reference_Allele, and Tumor_Seq_Allele1+2 columns
'''
# outfmaf = open( out+".maf", 'w')

flag_header = 1 
cnt_line = 0 
cnt_res = 0

with open(inp) as inpf: 
      for line in inpf.readlines():
            if re.match("^##",line):
		  outf.write(line) 
            elif re.match("^#C",line):
		  outf.write(line)
                  tempHeader = line.strip().split("\t")
                  sample1Name = tempHeader[-1]
            else:
                  cnt_line = cnt_line + 1
                  [chrom, pos, identifier, ref, alt, qual, filters, info, format,sample1] = \
                             line.strip().split("\t")
                  infoArray = info.split(";")

                  if not dbsnp in infoArray:
			infoArray_value = [a for a in infoArray if re.findall(r'=',a)]
                        # print str(cnt_line) + "\t line processing..." 
                        #print "\t".join(infoArray) + "\n"  
                        infoDict = {}
                        flag=1
                        for infosub in infoArray_value:
                             temp = infosub.split("=")
			     #print "\t".join(temp)
			     #[tag,val] = infosub.split("=")
                             if len(temp) == 2:
                                   infoDict[temp[0]] = temp[1]
                             elif len(temp) == 1:
                                   infoDict[temp[0]] = ''
			    # print temp[0] + ":" + str(temp[1])
			    # print len(infoDict.keys())
			##-----calling filtration
                        if ('MQ0' in infoDict.keys()) and (int(infoDict['MQ0']) > int(infoDict['DP']) * mq0_cut):
                             flag = 0
                             # print "fail MQ0"
                             
                        if ('DP' in infoDict.keys()) and int(infoDict['DP']) < dp_cut:
                             flag = 0
                             # print "fail DP"

                        if 'DP4' in infoDict.keys() : 
                          [dp3,dp4] = infoDict['DP4'].split(",")[2:] 
                          if (int(dp3) < dp3_cut) or (int(dp4) < dp4_cut ) or (int(dp3 + dp4) < dp34_cut) :
                             flag = 0
                             # print "fail DP4"
                             
                        if ('CLR' in infoDict.keys()) and int(infoDict['CLR']) < clr_cut :
                             flag = 0
                             # print "fail CLR"

                        #-----------------variant
			if infoDict.get('VariantType','-') == variantType :
			     flag = 0
                             
                        if infoDict.get('SomaticCategory','-') == snv_type :
                             flag = 0
                            # print "fail SomaticCategory"
                        if infoDict.get('functionalClass','-') == funcClass_type :
                             flag = 0
                            # print "fail SomaticCategory"
                        if infoDict.get('geneName','-') in fpgeneArray :
                             flag = 0
                             
                        formatArray = format.split(":")
                        #assume first one cancer, second one normal
                        sample1Array = sample1.split(":")
                        if 'AD' in formatArray: 
                             idxAD = formatArray.index('AD')
                             s1temp = sample1Array[idxAD].split(",")
                             if len(s1temp) == 0  :      
                               sample1Ref = '-'
                               sample1Alt = '-'
                             elif len(s1temp) == 1 : 
                               sample1Ref = sample1Array[idxAD].split(",")[0]
                               sample1Alt = '-'
                             else :
                               [sample1Ref, sample1Alt] = sample1Array[idxAD].split(",")[:2]
                        else :
                             flag = 0
                             # print "fail AD"
                        sample1PL_homoRef = sample1Array[formatArray.index('PL')].split(",")[0]  
                        if int(sample1PL_homoRef) < pltum_cut :
                             flag = 0
                             
                        ###-------------preparing for output 
                        # if flag_header == 1 :
                             # outf.write("chrom\tpos\tid\tref\talt\tqual\tfilter\tGeneName\tfunc\tfuncClass\tRef_" + \
                            # sample1Name +"\t" + "Alt_" + sample1Name + "\t" + "Ref_" + sample2Name +"\t" + "Alt_" + sample2Name + "\t" + \
                            # "1KG.score\tCOSMIC.FREQ\tESP5400.score\tCOSMIC.NonSynonymous.Freq\tCOSMIC.Synonymous.Freq\tljb_pp2.score\tDrugBank.INFO\n")
                             # flag_header = 0
                        if flag == 1:
                              # print "find..." +infoDict['geneName']
			      outf.write(line)
                              cnt_res = cnt_res + 1
                              # outf.write(chrom + "\t" +  pos + "\t"  + identifier + "\t" + ref + "\t" +  alt + "\t" +  qual + "\t" + filters + "\t" \
                              #         + infoDict['geneName'] + "\t"  +  infoDict.get('function','-') + "\t"+  infoDict.get('functionalClass','-') + "\t"  \
                              #         + sample1Ref + "\t" +  sample1Alt + "\t" +  sample2Ref + "\t" +  sample2Alt + "\t" \
                              #         + infoDict.get('1KG.score','-') + "\t" + infoDict.get('ESP5400.score','-') + "\t" \
                              #         + infoDict.get('COSMIC.FREQ','-') + "\t" +  infoDict.get('COSMIC.NonSynonymous.Freq','-') + "\t"\
                              #         + infoDict.get('COSMIC.Synonymous.Freq','-') + "\t" +  infoDict.get('ljb_pp2.score','-') + "\t" \
                              #         + infoDict.get('DrugBank.INFO','-')  \
                              #         +  "\n") 
                        else:
                              continue
                  else :
                       continue
      #infoDict['']
      #mafnorm
      #maftum
      #break
print "mutations processed: " + str(cnt_line) 
print "mutations survived after filtering: " + str(cnt_res)
print "#-----END------"

outf.close()

