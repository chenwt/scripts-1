#!/usr/bin/env python
#input: 1.<file:vcf after gatk,annovar,and mutable filtering process >
#	2.<file:potential false positive gene names, one by one line>
#	3.<string: output file name> 
#output: 1.<file: VCF after filtering> 
#desp: filtering each final vcfs 

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

usage = "Usage: filterIndel.py -i <input.vcf> -o <out.vcf>
	 Example: filterIndel.py -i /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/indelCall/test/TCGA-A1-A0SD-01A-11D-A10Y-09.bam_dindel_ouput.variantCalls.VCF -o /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/indelCall/test/TCGA-A1-A0SD-01A-11D-A10Y-09.bam_dindel_ouput.variantCalls.VCF.filtered.VCF"
argv = sys.argv[1:]
try:
      opts,args = getopt.getopt(argv,"hi:o:f:")
except getopt.GetoptError:
      print "Error Command line parameter \n" + usage 

for opt,arg in opts:
      if opt == "-h":
	    print usage 
      elif opt == "-i":
            inp = arg
      elif opt == "-o":
            out = arg
            outlog = out + ".log"

##1       1276973 .       G       GACAC   379     PASS    DP=32;NF=9;NR=0;NRS=11;NFS=0;HP=1       GT:GQ   1/1:30
###CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO=Depth of coverage; number of reads covering variants in the forward strand;number of reads covering....in the reverse strand; number of reads ... variants site in forward; number of reads ...variants site in reverse;    FORMAT  SAMPLE
#setting parameter
#gtnorm_cut = '0/0'
pltum_cut = 40
#plnorm_cut = 30
#clr_cut = 1
dp3_cut = 1
dp4_cut = 1
dp34_cut = 3
dp_cut = 10
#mq0_cut =  0.1
#dpnorm_cut = 1
#mafnorm_cut = 0.05
#maftum_cut = 0
snv_type = 'GERMLINE'
funcClass_type = 'synonymousSNV'
fpgeneArray = []
with open(fpgenefile) as fpgenef:
  for line in fpgenef.readlines():
    fpgeneArray.append(line.strip())
print "fpgene Num: " + str(len(fpgeneArray))

outf = open(out,'w')
#comment_re = re.compile(
flag_header = 1 
cnt_line = 0 
cnt_res = 0
with open(inp) as inpf: 
      for line in inpf.readlines():
            if re.match("^##",line):
                  pass
            elif re.match("^#C",line):
                  tempHeader = line.strip().split("\t")
                  sample1Name = tempHeader[-2]
                  sample2Name = tempHeader[-1]
            else:
                  cnt_line = cnt_line + 1
                  [chrom, pos, identifier, ref, alt, qual, filters, info, format,sample1,sample2] = \
                             line.strip().split("\t")
                  if filters.split(";")[-1] == 'PASS' :
                        # print str(cnt_line) + "\t line processing..." 
                        ### getting info
                        infoArray = info.split(";")
                        #print "\t".join(infoArray) + "\n"  
                        infoDict = {}
                        flag=1
                        for infosub in infoArray:
                             temp = infosub.split("=")
                       #print "\t".join(temp)
                        #[tag,val] = infosub.split("=")
                             if len(temp) == 2:
                                   infoDict[temp[0]] = temp[1]
                             elif len(temp) == 1:
                                   infoDict[temp[0]] = ''
                         # print temp[0] + ":" + str(temp[1])
                        # print len(infoDict.keys())
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
                        sample2Array = sample2.split(":")
                        if 'AD' in formatArray: 
                             idxAD = formatArray.index('AD')
                             s1temp = sample1Array[idxAD].split(",")
                             s2temp = sample2Array[idxAD].split(",")
                             if len(s1temp) == 0  :      
                               sample1Ref = '-'
                               sample1Alt = '-'
                             elif len(s1temp) == 0 : 
                               sample2Ref = '-'
                               sample2Alt = '-'
                             elif len(s1temp) == 1 : 
                               sample1Ref = sample1Array[idxAD].split(",")[0]
                               sample1Alt = '-'
                             elif len(s2temp) == 1:
                               sample2Ref = sample2Array[idxAD].split(",")[0]
                               sample2Alt = '-'
                             else :
                               [sample1Ref, sample1Alt] = sample1Array[idxAD].split(",")[:2]
                               [sample2Ref, sample2Alt] = sample2Array[idxAD].split(",")[:2]
                        else :
                             flag = 0
                             # print "fail AD"
                        #maftum = int(sample1Alt) / int(sample1Ref) 
                        if sample2Array[formatArray.index('GT')] != gtnorm_cut:
                             flag = 0 
                             # print "fail GT"
                        if int(sample2Ref) + int(sample2Alt) < dpnorm_cut:
                             flag = 0 
                             # print "fail dpnorm"
                             
                        if int(sample2Ref) > 0 :
                             mafnorm = int(sample2Alt) / int(sample2Ref) 
                             if mafnorm > mafnorm_cut:
                                    flag = 0 
                                    
                        sample1PL_homoRef = sample1Array[formatArray.index('PL')].split(",")[0]  
                        sample2PL_hete = sample2Array[formatArray.index('PL')].split(",")[1]  
                        if int(sample1PL_homoRef) < pltum_cut :
                             flag = 0
                             
                        if int(sample2PL_hete) < plnorm_cut: 
                             flag = 0 
                             
                        ###-------------preparing for output 
                        if flag_header == 1 :
                             outf.write("chrom\tpos\tid\tref\talt\tqual\tfilter\tGeneName\tfunc\tfuncClass\tRef_" + \
                            sample1Name +"\t" + "Alt_" + sample1Name + "\t" + "Ref_" + sample2Name +"\t" + "Alt_" + sample2Name + "\t" + \
                            "1KG.score\tCOSMIC.FREQ\tESP5400.score\tCOSMIC.NonSynonymous.Freq\tCOSMIC.Synonymous.Freq\tljb_pp2.score\tDrugBank.INFO\n")
                             flag_header = 0
                        if flag == 1:
                              # print "find..." +infoDict['geneName']
                              cnt_res = cnt_res + 1
                              outf.write(chrom + "\t" +  pos + "\t"  + identifier + "\t" + ref + "\t" +  alt + "\t" +  qual + "\t" + filters + "\t" \
                                      + infoDict['geneName'] + "\t"  +  infoDict.get('function','-') + "\t"+  infoDict.get('functionalClass','-') + "\t"  \
                                      + sample1Ref + "\t" +  sample1Alt + "\t" +  sample2Ref + "\t" +  sample2Alt + "\t" \
                                      + infoDict.get('1KG.score','-') + "\t" + infoDict.get('ESP5400.score','-') + "\t" \
                                      + infoDict.get('COSMIC.FREQ','-') + "\t" +  infoDict.get('COSMIC.NonSynonymous.Freq','-') + "\t"\
                                      + infoDict.get('COSMIC.Synonymous.Freq','-') + "\t" +  infoDict.get('ljb_pp2.score','-') + "\t" \
                                      + infoDict.get('DrugBank.INFO','-')  \
                                      +  "\n") 
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

