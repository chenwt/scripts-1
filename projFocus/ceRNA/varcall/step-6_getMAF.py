#!/usr/bin/python
#J.HE
#Desp.: get som matrix for given gene list from a list of fils  
#input: 1.<VCF files directory> 2.<grandularity: snp/interval/gene > 3.< types of overlap with genes> 4. < gene list: first colun gene> 
#output: <file: mat>  

import os,getopt
import re
import sys
import operator
from optparse import OptionParser
from subprocess import Popen,PIPE
from itertools import islice

usage = "python " + sys.argv[0] + " \
        -d <VCF file directory>  \
        -g <genelist, with exStart, exEnd, targeted> \
        -o <output .mat file> \
        -r <locType: gene/ tss>  \
        -k <keyType: maf/mut/>"
# vcfdir = '/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/mafs/test'
# output = 'test.out'
# genelist = '/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/geneSamples/cancerGenes_nature12912_novel.list.geneSingleStartEnd'
# keyType= 'maf' 
# locType = 'tss'
example = "python " + sys.argv[0] + " -d $cwd -o $output -l $locType -k $keyType -g $genelist "
MSG = "MESSAGE:\t"
ERR = "ERROR:\t"
SUC = "SUCESS\t" 

parser = OptionParser()
parser.add_option("-d", type="string", dest="vcfDir")
parser.add_option("-o", type="string", dest="output")
parser.add_option("-l", type="string", dest="locType")
parser.add_option("-g", type="string", dest="genelist")
parser.add_option("-k", type="string", dest="keyType")

options, args = parser.parse_args()
optionDict = vars(options)

vcfdir = optionDict['vcfDir'] 
output = optionDict['output'] 
genelist = optionDict['genelist'] 
keyType= optionDict['keyType'] 
locType = optionDict['locType'] 

class Mutation:
  ''' class store one record from VCF file, include: chrom, posStart, posEnd, gene1, gene2, DP, ADref, ADalt, '''
  def __init__(self, chrom, posStart, posEnd):
    self.chrom = chrom 
    self.posStart = long(posStart)
    self.posEnd  = long(posEnd)
    self.name = ":".join(map(str,[chrom,posStart,posEnd])) 
  def addVal(self, valueDict):
    self.vals = valueDict
  def dispMut(self,key):
    print "\t".join(map(str,[self.chrom,self.posStart, self.posEnd, self.vals[key]]) )
  def getMAF(self):
    key = 'DP4'
    valKey = self.vals.get(key,'NA').strip()
    tempDp4 = map(float,valKey.split(","))
    tempRef = tempDp4[1] + tempDp4[2] 
    tempAlt = tempDp4[-1] + tempDp4[-2] 
    if tempRef > 0:
        maf = tempAlt / (tempRef + tempAlt)
    elif tempRef == 0:
        maf = 1
    else:
        maf = 'NA'
    return maf        

class Gene:
  def __init__(self,name, chrom, posStart, posEnd, strand):
      self.name = name 
      self.chrom = chrom
      self.posStart = long(posStart)
      self.posEnd = long(posEnd)
      self.strand = strand
      self.tss = ":".join(map(str,[chrom,posStart,posEnd]))
  def hasMut(self, mutation):
      if self.name.lower() == mutation.vals.get('geneName','0').strip().lower():
          return 1
      else :
          return 0
  def nearMut(self, mutation):
      distCut = 1000000
      if self.chrom == mutation.chrom:
          if self.posStart > mutation.posEnd + distCut or self.posEnd < mutation.posStart - distCut :
              return 0
          else :
              return 1
      else :
          return 0
  def dispGene(self):
      print "\t".join(map(str,\
      [self.name,self.chrom,self.posStart, self.posEnd, self.strand]))
 
def newMutation(vcfRecord): 
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  H_LS-BH-A0E0-11A-13D-A128-09
    recArray = vcfRecord.strip().split("\t") 
    som = Mutation(recArray[0], int(recArray[1]), int(recArray[1]) + len(recArray[4]) )
    infoDict = dict(item.split("=",1) for item in  recArray[7].split(";") if "=" in item)
    formatDict = dict(zip(recArray[-2].split(":"),recArray[-1].split(":"))) 
    allDict = dict(infoDict, **formatDict)
    #print allDict
    som.addVal(allDict)
    return som
   
def sortedGeneList(geneDict):
    return(sorted(geneDict,key=operator.attrgetter('name')))

##----get all target genelist
geneDict = {}
with open(genelist) as f:
     for row in f.readlines() : 
         row = row.strip().split("\t")
         geneDict[Gene(row[0],row[1],row[2],row[3],row[4])] = ''
print "number of input genes: \t " + str(len(geneDict)) 

vcfs = [file for file in os.listdir(vcfdir) if re.findall(r"vcf$", file)]
cntFile = 1
keyType = keyType.lower()
locType = locType.lower()

os.chdir(vcfdir)

if keyType == 'mut' and locType == 'gene':  
    outputH = open(output,'w')
    for vcf in vcfs:
      print vcf
      geneDict = dict((k,0) for (k, v) in geneDict.items())
      with open(vcf, buffering = 10870) as f:
        for line in f.readlines():
	     if re.findall(r"^##", line):
		  continue
	     elif re.findall(r"^#", line):
		  tempBarcode = line.strip().split("\t")[-1]
	     else:
		  tempMut = newMutation(line)
		  for gene in geneDict.keys():
   	             if gene.hasMut(tempMut) :
   	                 geneDict[gene] = max(geneDict.get(gene,0),1)
   	                 continue ## one mutation only in one gene
   	             else :
   	                 continue
      geneList = sortedGeneList(geneDict) ##changed to list type              
      if cntFile == 1:
          gnames = [g.name for g in geneList]
          outputH.write("geneName\t" + "\t".join(sorted(gnames)) + "\n"  )
          valOut = tempBarcode
          for key in geneList:
               valOut = valOut + "\t" + str(geneDict[key])
          outputH.write(valOut + "\n" )
          outputH.flush()
      else :
          valOut = tempBarcode
          for key in geneList:
               valOut = valOut + "\t" + str(geneDict[key])
          outputH.write(valOut + "\n" )
          outputH.flush()
      cntFile = cntFile + 1   
    outputH.close() 	              
    print "num of vcf files:\t" + str(cntFile-1)  
    PYTHON = "~/tools/python/Python_current/python"
    SRCgetMutMat = "~/bin/trfile"
    output2 = output + ".mat" 
    pGenMat = Popen(  SRCgetMutMat +" " + output + " " + \
                    output2 , stdout=PIPE, shell=True)
    err = pGenMat.communicate()[1]
    if not err:  
        pRm = Popen("rm " + output, shell=True)
        pRm.communicate()
        print SUC 
    else:
        print ERR + "in Getting Matrix file.... "
	sys.exit()
elif keyType == 'maf' and locType == 'tss':
    outputH = open(output,'w')
    cntTemp = 0
    cntVCF = 0
    geneDictOld = geneDict
    for vcf in vcfs:
      print vcf
      cntVCF = cntVCF + 1
      with open(vcf, buffering = 100087000 ) as f:
        for line in f.readlines():
            geneDict = dict((k,[]) for (k, v) in geneDictOld.items())
            if re.findall(r"^##", line):
	            continue
            elif re.findall(r"^#", line):
	            tempBarcode = line.strip().split("\t")[-1]
            else:
                tempMut = newMutation(line)
                for gene in geneDict.keys():
                    if gene.nearMut(tempMut):
                        cntTemp = cntTemp + 1
                        if cntTemp % 10000 == 0 :
                            print str(cntTemp) +" pairs tested.."
                        geneDict[gene].append(tempMut)
                    else:
                        pass
            if cntFile == 1:
                  outputH.write("\t".join(["geneName","mutCode","MAF","barcode" ]) + "\n")
                  cntFile = cntFile + 1   
                  for g in geneDict.keys():
                       m = geneDict[g]
                       if len(m) > 0:
                            for mm in m :
                                temp_mutCode = ":".join(map(str,[mm.chrom,mm.posStart,mm.posEnd]))
                            valOut = "\t".join(map(str,[g.name,temp_mutCode,mm.getMAF(),tempBarcode])) 
                            outputH.write(valOut + "\n" )
                            outputH.flush()
            else :
                for g in geneDict.keys():
                    m = geneDict[g]
                    if len(m) > 0:
                        for mm in m :
                            temp_mutCode = ":".join(map(str,[mm.chrom,mm.posStart,mm.posEnd]))
                        valOut =  "\t".join(map(str,[g.name,temp_mutCode,mm.getMAF(),tempBarcode]))
                        outputH.write(valOut + "\n" )
                        outputH.flush()
    print "num of pairs:\t" + str(cntTemp)
    print "num of vcf files:\t" + str(cntVCF)  
    outputH.close() 	              
    sys.exit()
else:
    print MSG + locType + ":" + keyType +" combinations are Underdevelopment!"
    sys.exit()
print "##-----END---"

